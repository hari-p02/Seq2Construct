# C-VAE for Protein Sequence Modification: Technical Design Document

This document provides a comprehensive overview of the design, implementation, and training strategy for the Conditional Variational Autoencoder (C-VAE) developed to predict protein sequence modifications (specifically deletions).

## 1. Data Intelligence: Protein Embeddings

The foundation of the model is its ability to understand the "biological language" of protein sequences. Instead of using simple one-hot encoding, we utilize high-dimensional semantic embeddings.

### Embedding Model: Ankh-large
We use **Ankh-large**, a T5-based protein language model specifically trained on curated biological data.
*   **Dimensionality**: 1536 dimensions per sequence.
*   **Process**:
    1.  The raw amino acid sequence is tokenized.
    2.  The T5 encoder produces a "Last Hidden State" with a vector for every token in the sequence.
    3.  **Mean Pooling**: We take the average of all token vectors (excluding the EOS/padding tokens) to create a single 1536-dimensional "protein identity card" for the entire sequence.

### Design Choice: Why Ankh-large?
Standard LSTMs struggle to learn the complex folding and evolutionary properties of proteins from scratch. By using Ankh-large, we give the model a "pre-educated" brain that already understands which amino acids are chemically similar and which positions are likely important for stability.

---

## 2. Model Architecture: The C-VAE

The model is a **Conditional Variational Autoencoder**. Its goal is to take a protein embedding (the condition) and a modification pattern, and learn the underlying distribution that governs these modifications.

### Architecture Breakdown (Binary Setup)
We focus on 2 classes: **0 (Keep)** and **1 (Delete)**.

#### Encoder (The Pattern Recognizer)
*   **Input**: Protein Embedding ($X$) + Label Mask ($Y$)
*   **Structure**: 
    *   $X$ is processed via an MLP.
    *   $Y$ is processed via a Bidirectional LSTM to capture local sequence context.
*   **Output**: Two vectors: $\mu$ (mean) and $\sigma$ (standard deviation). These represent a "cloud" of probability in a 128-dimensional latent space.

#### Decoder (The Generative Artist)
*   **Input**: Sampled Latent Vector ($z$) + Protein Embedding ($X$)
*   **Structure**: 
    *   Fused $z$ and $X$ are passed through an MLP to initialize a **Generative LSTM**.
    *   The LSTM predicts tokens (0 or 1) sequentially.
*   **Output**: Logits for each position in the sequence.

---

## 3. The "Genius" in the Training Schedule

VAEs are notoriously difficult to train because the model often tries to "short-circuit" the math. We implemented three specific strategies to ensure high-quality learning.

### Strategy 1: KL Divergence & The 6-Epoch Warmup
**What is KL Divergence?** It’s a penalty that forces the latent space to look like a standard Bell Curve ($N(0,1)$). Without it, the model would just memorize every protein perfectly (which is useless for generating new ones).

*   **Warmup (Epochs 1-6)**: We set the `KL_weight` to `0.0`. 
*   **Rationale**: We want the model to first master the "Copying" task. By the time we start the "Variational" training in Epoch 7, the Decoder already knows how to translate a latent vector into a sequence flawlessly. This prevents the "Posterior Collapse" where the model ignores the latent space entirely.

### Strategy 2: Scheduled Sampling (Teacher Forcing)
**What is Teacher Forcing?** During training, we usually give the model the "Answer Key" for token $t-1$ when predicting token $t$. 

*   **The Problem**: In the real world, the model won't have the answer key. If it only trains with the key, it will crash and burn the moment it makes a single mistake.
*   **The Schedule**: We start at `tf_ratio = 1.0` (all answers) and decay down to `tf_ratio = 0.1` (only 10% answers).
*   **Rationale**: This slowly "pulls the training wheels off." It forces the model to learn how to recover from its own small errors, leading to a much more robust model for inference.

### Strategy 3: Class Weighting
*   **The Imbalance**: In your data, "Keep" tokens (0) are much more common than "Delete" tokens (1).
*   **The Fix**: We calculate the inverse frequency and assign a higher "cost" to missing a Deletion.
*   **Rationale**: This ensures the model doesn't just get "lazy" and predict 0 for everything to get easy points.

---

## 4. Performance Engineering (The Speedup)

To leverage the **NVIDIA A100**, we implemented three "Turbo" features:
1.  **Batch Size (256)**: Massively increases throughput by parallelizing the math.
2.  **Mixed Precision (FP16)**: Uses 16-bit math for speed without losing accuracy.
3.  **torch.compile**: Performs "Kernel Fusion," merging multiple LSTM operations into a single fast CUDA execution. This reduced epoch time from ~1.5 hours to **~20 minutes**.

---

## 5. Conclusion: Why 2 Classes (0/1)?

By merging "Class 2" (Modify) into "Class 0" (Keep), we eliminated the extreme **45x imbalance** that was causing the model to hallucinate modifications. The result is a model that is significantly more stable, has much lower loss, and provides cleaner deletion maps for your protein constructs.
