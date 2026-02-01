#!/usr/bin/env python3
"""
Protein Construct Prediction - Inference Script
Generate construct suggestions from protein sequences using trained C-VAE model
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
import sys
from pathlib import Path
from typing import List, Tuple
import argparse

# Add color support for terminal output
class Colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


# ============================================
# Model Architecture (must match training)
# ============================================

class Encoder(nn.Module):
    """Encoder: processes embeddings and masks to produce latent distribution"""
    
    def __init__(self, embedding_dim=1536, hidden_dim=512, lstm_hidden=256, latent_dim=128, dropout=0.2):
        super().__init__()
        
        # Embedding branch (MLP)
        self.embedding_mlp = nn.Sequential(
            nn.Linear(embedding_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.LayerNorm(hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout)
        )
        
        # Mask branch (LSTM)
        self.mask_embedding = nn.Embedding(3, 32)
        self.mask_lstm = nn.LSTM(
            input_size=32,
            hidden_size=lstm_hidden,
            num_layers=2,
            batch_first=True,
            dropout=dropout,
            bidirectional=True
        )
        
        # Fusion layer
        fusion_input_dim = hidden_dim // 2 + lstm_hidden * 2
        self.fusion = nn.Sequential(
            nn.Linear(fusion_input_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout)
        )
        
        # Output layers
        self.fc_mu = nn.Linear(hidden_dim, latent_dim)
        self.fc_logvar = nn.Linear(hidden_dim, latent_dim)
        
    def forward(self, embeddings, masks, attention_mask):
        emb_features = self.embedding_mlp(embeddings)
        mask_embedded = self.mask_embedding(masks)
        
        lengths = attention_mask.sum(dim=1).cpu()
        packed_input = nn.utils.rnn.pack_padded_sequence(
            mask_embedded, lengths, batch_first=True, enforce_sorted=False
        )
        
        packed_output, (hidden, cell) = self.mask_lstm(packed_input)
        
        hidden_fwd = hidden[-2]
        hidden_bwd = hidden[-1]
        mask_features = torch.cat([hidden_fwd, hidden_bwd], dim=1)
        
        combined = torch.cat([emb_features, mask_features], dim=1)
        fused = self.fusion(combined)
        
        mu = self.fc_mu(fused)
        logvar = self.fc_logvar(fused)
        
        return mu, logvar


class Decoder(nn.Module):
    """Decoder: generates modification mask sequence from latent vector and embedding"""
    
    def __init__(self, embedding_dim=1536, hidden_dim=512, lstm_hidden=256, latent_dim=128, num_classes=3, dropout=0.2):
        super().__init__()
        
        self.lstm_hidden = lstm_hidden
        self.num_classes = num_classes
        
        # Embedding branch (MLP)
        self.embedding_mlp = nn.Sequential(
            nn.Linear(embedding_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.LayerNorm(hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout)
        )
        
        # Fusion with latent vector
        fusion_input_dim = hidden_dim // 2 + latent_dim
        self.fusion = nn.Sequential(
            nn.Linear(fusion_input_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout)
        )
        
        # Project to LSTM initial hidden state
        self.fc_hidden = nn.Linear(hidden_dim, lstm_hidden)
        self.fc_cell = nn.Linear(hidden_dim, lstm_hidden)
        
        # Generative LSTM
        self.input_embedding = nn.Embedding(num_classes, 32)
        self.lstm = nn.LSTM(
            input_size=32 + hidden_dim,
            hidden_size=lstm_hidden,
            num_layers=2,
            batch_first=True,
            dropout=dropout
        )
        
        # Output layer
        self.output_layer = nn.Linear(lstm_hidden, num_classes)
        
    def forward(self, embeddings, z, target_masks, attention_mask, teacher_forcing_ratio=0.5):
        batch_size, seq_len = target_masks.shape
        
        emb_features = self.embedding_mlp(embeddings)
        combined = torch.cat([emb_features, z], dim=1)
        context = self.fusion(combined)
        
        h0 = self.fc_hidden(context).unsqueeze(0).repeat(2, 1, 1)
        c0 = self.fc_cell(context).unsqueeze(0).repeat(2, 1, 1)
        
        outputs = []
        hidden = (h0, c0)
        input_token = torch.zeros(batch_size, dtype=torch.long, device=embeddings.device)
        
        for t in range(seq_len):
            input_embedded = self.input_embedding(input_token)
            lstm_input = torch.cat([input_embedded, context], dim=1).unsqueeze(1)
            lstm_out, hidden = self.lstm(lstm_input, hidden)
            logits = self.output_layer(lstm_out.squeeze(1))
            outputs.append(logits)
            
            if t < seq_len - 1:
                if np.random.random() < teacher_forcing_ratio:
                    input_token = target_masks[:, t]
                else:
                    input_token = logits.argmax(dim=1)
        
        outputs = torch.stack(outputs, dim=1)
        return outputs


class ConditionalVAE(nn.Module):
    """Complete Conditional VAE model"""
    
    def __init__(self, embedding_dim=1536, hidden_dim=512, lstm_hidden=256, 
                 latent_dim=128, num_classes=3, dropout=0.2):
        super().__init__()
        
        self.encoder = Encoder(embedding_dim, hidden_dim, lstm_hidden, latent_dim, dropout)
        self.decoder = Decoder(embedding_dim, hidden_dim, lstm_hidden, latent_dim, num_classes, dropout)
        
    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std
    
    def forward(self, embeddings, masks, attention_mask, teacher_forcing_ratio=0.5):
        mu, logvar = self.encoder(embeddings, masks, attention_mask)
        z = self.reparameterize(mu, logvar)
        logits = self.decoder(embeddings, z, masks, attention_mask, teacher_forcing_ratio)
        return logits, mu, logvar


# ============================================
# Embedding Generation (Ankh-large)
# ============================================

# Global variables for embedding model (loaded once)
_embedding_model = None
_embedding_tokenizer = None
_embedding_device = None

def load_embedding_model(device: torch.device):
    """Load Ankh-large model for generating embeddings (called once)"""
    global _embedding_model, _embedding_tokenizer, _embedding_device
    
    if _embedding_model is not None:
        return  # Already loaded
    
    print(f"{Colors.CYAN}Loading Ankh-large embedding model...{Colors.END}")
    
    from transformers import T5EncoderModel, AutoTokenizer
    
    hf_token = ""
    model_name = "ElnaggarLab/ankh-large"
    
    try:
        _embedding_tokenizer = AutoTokenizer.from_pretrained(
            model_name,
            do_lower_case=False,
            token=hf_token
        )
        
        _embedding_model = T5EncoderModel.from_pretrained(
            model_name,
            token=hf_token
        )
        
        _embedding_model.to(device).eval()
        _embedding_device = device
        
        print(f"{Colors.GREEN}✓ Ankh-large model loaded successfully!{Colors.END}\n")
        
    except Exception as e:
        print(f"{Colors.RED}❌ Error loading Ankh-large model: {e}{Colors.END}")
        print(f"{Colors.YELLOW}Falling back to placeholder embeddings...{Colors.END}\n")
        _embedding_model = None


def generate_embedding(sequence: str) -> np.ndarray:
    """
    Generate 1536-dim embedding for a protein sequence using Ankh-large.
    
    Uses mean pooling over the sequence (excluding EOS token).
    """
    global _embedding_model, _embedding_tokenizer, _embedding_device
    
    # Fallback to random if model not loaded
    if _embedding_model is None:
        np.random.seed(hash(sequence) % (2**32))
        return np.random.randn(1536).astype(np.float32)
    
    try:
        # Tokenize
        inputs = _embedding_tokenizer(
            sequence,
            truncation=True,
            max_length=2048,
            padding=False,
            return_tensors="pt"
        )
        
        inputs = {k: v.to(_embedding_device) for k, v in inputs.items()}
        
        # Generate embeddings
        with torch.no_grad():
            outputs = _embedding_model(**inputs)
        
        # Mean pooling (exclude EOS token at the end)
        embeddings = outputs.last_hidden_state
        protein_vector = torch.mean(embeddings[0, :-1, :], dim=0)
        
        # Convert to numpy
        embedding = protein_vector.cpu().numpy().astype(np.float32)
        
        return embedding
        
    except Exception as e:
        print(f"{Colors.RED}❌ Error generating embedding: {e}{Colors.END}")
        print(f"{Colors.YELLOW}Using random fallback...{Colors.END}")
        np.random.seed(hash(sequence) % (2**32))
        return np.random.randn(1536).astype(np.float32)


# ============================================
# Inference Functions
# ============================================

def load_model(checkpoint_path: str, device: torch.device) -> ConditionalVAE:
    """Load trained model from checkpoint"""
    
    # Initialize model with architecture matching training
    model = ConditionalVAE(
        embedding_dim=1536,
        hidden_dim=512,
        lstm_hidden=256,
        latent_dim=128,
        num_classes=3,
        dropout=0.2
    )
    
    # Load checkpoint
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
    
    # Handle different checkpoint formats
    if 'model_state_dict' in checkpoint:
        state_dict = checkpoint['model_state_dict']
        # Remove torch.compile prefix if present
        if any(k.startswith('_orig_mod.') for k in state_dict.keys()):
            state_dict = {k.replace('_orig_mod.', ''): v for k, v in state_dict.items()}
        model.load_state_dict(state_dict)
    else:
        model.load_state_dict(checkpoint)
    
    model.to(device)
    model.eval()
    
    return model


def predict_constructs(
    model: ConditionalVAE,
    sequence: str,
    device: torch.device,
    num_samples: int = 5,
    temperature: float = 1.0,
    logit_bias: List[float] = None
) -> List[Tuple[str, np.ndarray, float]]:
    """
    Generate construct suggestions for a protein sequence.
    
    Args:
        logit_bias: List of 3 floats to add to [Keep, Delete, Modify] logits
                   to calibrate sensitivity.
    """
    
    # Generate embedding
    embedding = generate_embedding(sequence)
    embedding_tensor = torch.tensor(embedding, dtype=torch.float32).unsqueeze(0).to(device)
    
    # Create dummy input mask (all zeros - no modifications)
    seq_len = len(sequence)
    dummy_mask = torch.zeros(1, seq_len, dtype=torch.long).to(device)
    attention_mask = torch.ones(1, seq_len, dtype=torch.bool).to(device)
    
    constructs = []
    
    with torch.no_grad():
        for i in range(num_samples):
            # Forward pass (no teacher forcing)
            logits, mu, logvar = model(embedding_tensor, dummy_mask, attention_mask, teacher_forcing_ratio=0.0)
            
            # Apply logit bias to handle class imbalance if provided
            if logit_bias is not None:
                bias_tensor = torch.tensor(logit_bias, device=device).view(1, 1, 3)
                logits = logits + bias_tensor
                
            # Apply temperature scaling
            logits = logits / temperature
            
            # Get predictions
            probs = F.softmax(logits, dim=-1)
            predictions = logits.argmax(dim=-1).squeeze(0).cpu().numpy()
            
            # Calculate confidence (average max probability)
            confidence = probs.max(dim=-1)[0].mean().item()
            
            # Generate construct sequence
            construct = generate_construct(sequence, predictions)
            
            constructs.append((construct, predictions, confidence))
    
    # Sort by confidence
    constructs.sort(key=lambda x: x[2], reverse=True)
    
    return constructs


def generate_construct(sequence: str, modification_mask: np.ndarray) -> str:
    """
    Generate construct sequence from original sequence and modification mask.
    
    Mask values:
        0: Keep amino acid
        1: Delete amino acid
        2: Modify amino acid (placeholder for now)
    """
    construct = []
    for i, (aa, mod) in enumerate(zip(sequence, modification_mask)):
        if mod == 0:
            construct.append(aa)
        elif mod == 1:
            # Skip (deletion)
            continue
        elif mod == 2:
            # Modification (for now, mark with lowercase)
            construct.append(aa.lower())
    
    return ''.join(construct)


# ============================================
# Visualization
# ============================================

def print_header():
    """Print fancy header"""
    print(f"\n{Colors.CYAN}{Colors.BOLD}{'='*70}{Colors.END}")
    print(f"{Colors.CYAN}{Colors.BOLD}  🧬 Protein Construct Prediction - C-VAE Inference{Colors.END}")
    print(f"{Colors.CYAN}{Colors.BOLD}{'='*70}{Colors.END}\n")


def print_construct(idx: int, construct: str, mask: np.ndarray, confidence: float, original: str):
    """Print a single construct suggestion with visual formatting"""
    
    print(f"{Colors.BOLD}{Colors.GREEN}Construct #{idx}{Colors.END}")
    print(f"{Colors.BOLD}Confidence:{Colors.END} {confidence*100:.2f}%")
    print(f"{Colors.BOLD}Length:{Colors.END} {len(construct)} aa (original: {len(original)} aa)")
    
    # Count modifications
    deletions = np.sum(mask == 1)
    modifications = np.sum(mask == 2)
    kept = np.sum(mask == 0)
    
    print(f"{Colors.BOLD}Modifications:{Colors.END} {kept} kept, {deletions} deleted, {modifications} modified")
    
    # Print sequence with color coding
    print(f"\n{Colors.BOLD}Original:{Colors.END}")
    print_colored_sequence(original, mask)
    
    print(f"\n{Colors.BOLD}Construct:{Colors.END}")
    print(f"{Colors.GREEN}{construct}{Colors.END}")
    
    print(f"\n{Colors.CYAN}{'─'*70}{Colors.END}\n")


def print_colored_sequence(sequence: str, mask: np.ndarray):
    """Print sequence with color-coded modifications"""
    output = []
    for aa, mod in zip(sequence, mask):
        if mod == 0:
            output.append(f"{Colors.GREEN}{aa}{Colors.END}")  # Keep - green
        elif mod == 1:
            output.append(f"{Colors.RED}{aa}{Colors.END}")    # Delete - red
        elif mod == 2:
            output.append(f"{Colors.YELLOW}{aa}{Colors.END}") # Modify - yellow
    
    # Print in chunks of 60
    seq_str = ''.join(output)
    for i in range(0, len(output), 60):
        chunk = ''.join(output[i:i+60])
        print(f"  {chunk}")


# ============================================
# Main Interface
# ============================================

def interactive_mode(model: ConditionalVAE, device: torch.device):
    """Interactive terminal interface"""
    
    print_header()
    print(f"{Colors.BOLD}Interactive Mode{Colors.END}")
    print(f"Enter protein sequences to get construct suggestions.")
    print(f"Type 'quit' or 'exit' to stop.\n")
    
    while True:
        # Get input
        print(f"{Colors.BOLD}Enter protein sequence:{Colors.END} ", end='')
        sequence = input().strip().upper()
        
        if sequence.lower() in ['quit', 'exit', 'q']:
            print(f"\n{Colors.CYAN}Goodbye! 👋{Colors.END}\n")
            break
        
        if not sequence:
            continue
        
        # Validate sequence
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        if not all(aa in valid_aa for aa in sequence):
            print(f"{Colors.RED}❌ Invalid sequence! Use only standard amino acids.{Colors.END}\n")
            continue
        
        print(f"\n{Colors.CYAN}Generating constructs...{Colors.END}\n")
        
        # Generate predictions
        constructs = predict_constructs(model, sequence, device, num_samples=5)
        
        # Display results
        for i, (construct, mask, confidence) in enumerate(constructs, 1):
            print_construct(i, construct, mask, confidence, sequence)
        
        print()


def main():
    parser = argparse.ArgumentParser(
        description='Generate protein construct suggestions using trained C-VAE model',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Interactive mode
  python inference.py
  
  # Single sequence
  python inference.py --sequence "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL"
  
  # Specify number of suggestions
  python inference.py --sequence "MKTAYIAK..." --num-samples 10
        """
    )
    
    parser.add_argument(
        '--model', '-m',
        type=str,
        default='/Users/haripat/Desktop/SF/protein/models/v1.pt',
        help='Path to model checkpoint (default: v1.pt)'
    )
    
    parser.add_argument(
        '--sequence', '-s',
        type=str,
        help='Protein sequence to analyze (if not provided, enters interactive mode)'
    )
    
    parser.add_argument(
        '--num-samples', '-n',
        type=int,
        default=5,
        help='Number of construct suggestions to generate (default: 5)'
    )
    
    parser.add_argument(
        '--temperature', '-t',
        type=float,
        default=1.0,
        help='Sampling temperature (higher = more diverse, default: 1.0)'
    )
    
    parser.add_argument(
        '--logit-bias',
        type=float,
        nargs=3,
        default=None,
        help='Bias for [Keep, Delete, Modify] logits (e.g., 2.0 0.0 0.0 to favor Keeping)'
    )
    
    parser.add_argument(
        '--keep-bias',
        type=float,
        default=0.0,
        help='Shortcut to add bias specifically to the "Keep" class (default: 0.0)'
    )
    
    args = parser.parse_args()
    
    # Process biases
    logit_bias = args.logit_bias
    if logit_bias is None and args.keep_bias != 0.0:
        logit_bias = [args.keep_bias, 0.0, 0.0]
    
    # Setup device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    # Print header
    print_header()
    print(f"{Colors.BOLD}Device:{Colors.END} {device}")
    print(f"{Colors.BOLD}Model:{Colors.END} {args.model}")
    if logit_bias:
        print(f"{Colors.BOLD}Logit Bias:{Colors.END} {logit_bias}")
    
    # Load model
    try:
        print(f"\n{Colors.CYAN}Loading model...{Colors.END}")
        model = load_model(args.model, device)
        print(f"{Colors.GREEN}✓ Model loaded successfully!{Colors.END}\n")
    except Exception as e:
        print(f"{Colors.RED}❌ Error loading model: {e}{Colors.END}")
        sys.exit(1)
    
    # Run inference
    if args.sequence:
        # Single sequence mode
        sequence = args.sequence.strip().upper()
        
        print(f"{Colors.BOLD}Input sequence:{Colors.END} {sequence[:60]}{'...' if len(sequence) > 60 else ''}")
        print(f"{Colors.BOLD}Length:{Colors.END} {len(sequence)} amino acids\n")
        
        print(f"{Colors.CYAN}Generating {args.num_samples} construct suggestions...{Colors.END}\n")
        
        constructs = predict_constructs(
            model, sequence, device, 
            num_samples=args.num_samples, 
            temperature=args.temperature,
            logit_bias=logit_bias
        )
        
        for i, (construct, mask, confidence) in enumerate(constructs, 1):
            print_construct(i, construct, mask, confidence, sequence)
    else:
        # Interactive mode
        interactive_mode(model, device)


if __name__ == '__main__':
    main()
