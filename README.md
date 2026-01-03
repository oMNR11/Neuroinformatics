# Neuroinformatics
This is a repository for uploading documents and codes (primarily in MATLAB) for the course contents of Neuroinformatics.

I will be following Mike X Cohen's Analyzing Neural Time Series Data
```mermaid
graph TD
    %% 1. Raw Data & Preprocessing
    subgraph Preprocessing ["Step 1: Signal Processing"]
        Raw[("Raw ECoG Data\n(Continuous)")] --> Filter["Bandpass Filter\n(5 Bands: δ, θ, α, β, γ)"]
        Filter --> Hilbert["Hilbert Transform"]
        Hilbert --> Power[("Continuous Oscillatory Power\n(All Electrodes)")]
    end

    %% 2. Epoching
    subgraph Epoching ["Step 2: Segmentation (Epoching)"]
        Power --> CutStudy["Cut STUDY Events\n(200 to 1600ms post-stimulus)"]
        Power --> CutRecall["Cut RECALL Events\n(-600 to 200ms relative to voice)"]
        CutStudy --> Stack["Stack into Single Matrix\n(Rows = Events, Cols = Elec x Freq)"]
        CutRecall --> Stack
    end

    %% 3. Dimensionality Reduction
    subgraph DimRed ["Step 3: Dimensionality Reduction"]
        Stack --> ZScore["Z-Transform Power"]
        ZScore --> PCA["Run PCA\n(Principal Component Analysis)"]
        PCA --> PCs[("Principal Components\n(Neural 'States')")]
    end

    %% 4. Feature Selection
    subgraph Selection ["Step 4: Context Feature Selection"]
        PCs --> ExtractStudy["Isolate STUDY Rows only"]
        ExtractStudy --> AutoCorr{"Calculate Lag-1\nAutocorrelation (ρ)"}
        AutoCorr -- "ρ > 0 (Slow Drift)" --> Keep["Keep Component\n(Context Feature)"]
        AutoCorr -- "ρ ≤ 0 (Noise/Transient)" --> Discard["Discard Component"]
    end

    %% 5. Similarity Analysis
    subgraph Analysis ["Step 5: Lag Analysis"]
        Keep --> VecStudy["Vector f_study\n(Target Pattern)"]
        Keep --> VecRecall["Vector f_recall\n(Probe Pattern)"]
        
        VecRecall -- "Compare vs." --> DotProd("Normalized Dot Product\n(Similarity Calculation)")
        VecStudy --> DotProd
        
        DotProd --> LagAssign["Assign Lag\n(0 = Same Item, ±1 = Neighbors)"]
    end

    %% 6. Final Result
    subgraph Result ["Step 6: Output"]
        LagAssign --> Binning["Average Similarity\nper Lag Bin"]
        Binning --> Graph[("Final Graph:\nNeural Similarity vs. Lag")]
    end
```
