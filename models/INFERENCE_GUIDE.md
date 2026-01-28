# Protein Construct Prediction - Inference Guide

## Quick Start

### Interactive Mode (Recommended)
```bash
cd /Users/haripat/Desktop/SF/protein/models
python3 inference.py
```

Then enter protein sequences when prompted. Type `quit` to exit.

### Single Sequence Mode
```bash
python3 inference.py --sequence "MKTAYIAKQRQISFVKSHFSRQLE"
```

### Advanced Options
```bash
# Generate 10 suggestions
python3 inference.py --sequence "MKTAYIAK..." --num-samples 10

# Use different model checkpoint
python3 inference.py --model best_model.pt --sequence "MKTAYIAK..."

# Adjust sampling temperature (higher = more diverse)
python3 inference.py --sequence "MKTAYIAK..." --temperature 1.5
```

## Output Explanation

The script generates **construct suggestions** with:

- **Confidence Score**: Model's confidence in the prediction (0-100%)
- **Modifications**: 
  - 🟢 **Green**: Amino acids kept in construct
  - 🔴 **Red**: Amino acids deleted
  - 🟡 **Yellow**: Amino acids modified
- **Construct Sequence**: Final suggested protein construct

## Example Output

```
🧬 Protein Construct Prediction - C-VAE Inference
======================================================================

Construct #1
Confidence: 87.34%
Length: 180 aa (original: 244 aa)
Modifications: 180 kept, 64 deleted, 0 modified

Original:
  MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVK...

Construct:
  MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVK...
```

## Important Notes

⚠️ **Embedding Model Required**: The current script uses **placeholder embeddings**. For production use, you need to:

1. Install your embedding model (e.g., ESM-2, ProtBERT)
2. Replace the `generate_embedding()` function in `inference.py`

Example with ESM-2:
```python
from transformers import AutoTokenizer, AutoModel

tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D")
model = AutoModel.from_pretrained("facebook/esm2_t33_650M_UR50D")

def generate_embedding(sequence: str) -> np.ndarray:
    inputs = tokenizer(sequence, return_tensors="pt")
    with torch.no_grad():
        outputs = model(**inputs)
        # Mean pooling
        embedding = outputs.last_hidden_state.mean(dim=1).squeeze().numpy()
    return embedding
```

## Modification Mask Interpretation

- **0**: Keep amino acid (no change)
- **1**: Delete amino acid (remove from construct)
- **2**: Modify amino acid (currently shown as lowercase)

## Tips

- Longer sequences may take more time to process
- Higher `--num-samples` gives more diverse suggestions
- Higher `--temperature` increases diversity but may reduce quality
- The model works best on sequences similar to training data

## Troubleshooting

**Error: "CUDA device not available"**
- The script automatically falls back to CPU
- No action needed

**Error: "Invalid sequence"**
- Use only standard 20 amino acids (ACDEFGHIKLMNPQRSTVWY)
- Remove any special characters or numbers

**Low confidence scores**
- Try adjusting temperature
- Sequence may be very different from training data
- Consider retraining with more similar examples
