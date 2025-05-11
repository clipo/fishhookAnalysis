# Fishhook Morphometric Analysis

## Overview

This Shiny application provides a comprehensive platform for the morphometric analysis of archaeological or modern fishhooks based on digital images. The application implements a series of advanced computational techniques for outline extraction, shape analysis, and functional region identification to facilitate quantitative comparative studies of fishhook morphology.

## Features

### Image Processing
- **Multiple Image Sources**: Use included demo images or upload your own JPEG photographs
- **Interactive Cropping**: Select regions of interest within images
- **Adaptive Thresholding**: Configure threshold levels for artifacts that are darker or lighter than background
- **Noise Reduction**: Remove small artifacts and improve outline quality

### Fishhook Region Detection
Three different methods for identifying functional regions:

1. **Geometric Method**: Uses curvature analysis to automatically identify key regions based on shape features
2. **Proportional Method**: Defines regions based on user-specified percentages of outline points
3. **Manual Method**: Interactive selection of key points on the outline

### Morphometric Analysis
- **Outline Extraction**: Automatic contour detection and interpolation
- **Elliptical Fourier Analysis (EFA)**: Configurable number of harmonics
- **Principal Component Analysis (PCA)**: Visualization of shape variation in morphospace
- **Clustering**: Automatic grouping of similar shapes

### Functional Metrics
Automatically calculates critical functional dimensions:
- **Gape Width**: Distance between shank and point
- **Throat Depth**: Maximum perpendicular distance from bend to shank-point line
- **Point Angle**: Angular measurement at the point tip

### Visualization
- **Grid View**: Display all processed outlines
- **Stacked View**: Overlaid aligned outlines for direct comparison
- **Harmonic Power**: Analysis of Fourier harmonic contribution
- **PCA Morphospace**: 2D plot with visualization of shape extremes
- **Thin-Plate Splines**: Visualization of shape deformation along PC axes
- **Regional Analysis**: Color-coded visualization of fishhook functional regions

### Data Export
- **PDF Report**: Generate comprehensive analysis report with visualizations

## Installation

### Prerequisites
- R (version 4.0 or higher recommended)
- Required R packages:
  - shiny
  - Momocs
  - imager
  - dplyr
  - ggplot2

### Installation Steps

1. Install required R packages:
```R
install.packages(c("shiny", "Momocs", "imager", "dplyr", "ggplot2"))
```

2. Clone or download this repository:
```bash
git clone https://github.com/clipo/fishhookAnalysis.git
cd fishhookAnalysis
```

3. Run the application:
```R
shiny::runApp()
```

## Getting Started: A Detailed Guide

### Image Requirements

For optimal results, prepare your fishhook images as follows:

- **Format**: JPEG format (.jpg or .jpeg)
- **Resolution**: At least 300 DPI recommended (higher is better)
- **Background**: High contrast with the fishhook (white or black background works best)
- **Orientation**: Fishhook should be oriented with the point facing right or down
- **Lighting**: Even lighting without strong shadows or reflections
- **Scale**: Include a scale bar in images if absolute size is important

### Demo Images

To get started quickly:
1. Place sample fishhook JPEG images in the 'www' folder
2. Select "Demo Images" in the application interface

### Step-by-Step Workflow

#### 1. Prepare and Load Images
- **For demo images**: Place JPEG files in the 'www' folder and select "Demo Images"
- **For your own images**: Select "Upload My Own" and use the file upload dialog

#### 2. Process Each Image
- The first image is displayed automatically
- **Image Cropping**:
  - Click two points (opposite corners) on the image to define a crop rectangle
  - Crop should include the entire fishhook with minimal background
- **Configure Settings**:
  - Set "Fishhook appearance" to "Darker than background" or "Lighter than background"
  - Adjust "Threshold" slider until the preview shows a clean binary representation
  - Use "Noise removal" to eliminate small artifacts
- **Process the Image**:
  - Click "Process Image" to extract the outline
  - The outlined fishhook appears in the bottom panel
- **Move to Next Image**:
  - Click "Next Image" to proceed to the next fishhook
  - Repeat the process for all images

#### 3. Configure Analysis Settings
- **Harmonics**: Set the number of harmonics for EFA (typically 20-30)
- **Clusters**: Set the number of clusters for shape classification

#### 4. Define Regions
Choose one of three region detection methods:

- **Geometric (Shape-based)**:
  - Uses curvature analysis to automatically identify regions
  - Adjust "Curvature sensitivity" to control detection
  - Adjust "Curvature smoothing" to reduce noise
- **Proportional (Percentage-based)**:
  - Manually set the percentage of points for each region
  - Shank, bend, point, and barb proportions can be adjusted
- **Manual key points**:
  - Process an image first
  - Click buttons to define key points on the outline:
    1. "Top of Shank" - upper end of shank
    2. "Bottom of Shank" - lower end of shank
    3. "Deepest Bend Point" - deepest part of bend
    4. "Point Tip" - tip of the hook point
    5. "Barb" - barb tip (if present)

#### 5. Run Analysis
- Click "Run Analysis" to perform calculations
- The application will:
  - Align and normalize all outlines
  - Perform Elliptical Fourier Analysis
  - Calculate principal components
  - Identify clusters
  - Calculate functional metrics
  - Generate visualizations

#### 6. Interpret Results

The "Results" tab contains multiple visualizations:

- **Outlines Grid**: Shows all processed outlines with labels
- **Stacked Outlines**: All outlines superimposed and aligned
- **EFA Harmonic Power**: 
  - Shows the contribution of each harmonic
  - The bars represent each harmonic's power percentage
  - The red line shows cumulative power
  - Higher percentage for early harmonics indicates simpler shapes
- **PCA Morphospace**:
  - Points represent individual fishhooks
  - Colors indicate assigned clusters
  - Axes represent principal components of shape variation
  - PC1 captures the most shape variation, PC2 the second most
- **Thin-Plate Splines**:
  - Shows shape deformation along PC axes
  - Visualizes how shapes differ at extremes of PCA
- **Integration/Modularity**:
  - Color-coded regions of the fishhook
  - Functional metrics visualization
  - Shows relationships between different parts

#### 7. Fishhook Metrics

The "Fishhook Metrics" panel displays:
- **Gape width**: Distance between shank and point
- **Throat depth**: Maximum depth of bend region
- **Point angle**: Angle at the tip of the point

#### 8. Export Results
- Click "Download Analysis (PDF)" to generate a comprehensive report
- The PDF includes:
  - Summary of analysis parameters
  - All visualizations and results
  - Metrics for analyzed fishhooks

## Practical Tips

### Optimizing Image Processing
- Start with a threshold of 0.5 and adjust up or down based on preview
- Increase noise removal if small artifacts appear in the binary image
- If outline extraction fails, try adjusting lighting/contrast in the original image

### Region Detection Methods
- **Geometric method** works best with clear, well-defined shapes
- **Proportional method** is more reliable for irregular or eroded fishhooks
- **Manual method** offers the most control but requires more time

### Batch Processing
- Process all images first before running analysis
- Use "Next Image" to move through your dataset
- The "Reset" button clears current processing if you need to start over

### Common Issues and Solutions

| Issue | Solution |
|-------|----------|
| No outline detected | Adjust threshold, increase contrast in original image |
| Noisy outline | Increase noise removal slider |
| Missing parts of hook | Adjust cropping to include entire hook |
| Incorrect region detection | Try a different detection method |
| Analysis errors | Increase number of points using manual method |

## Interpreting Outputs

### Shape Similarity
- Hooks that cluster together in PCA morphospace have similar shapes
- Distance between points corresponds to degree of shape difference

### Functional Interpretation
- Gape width relates to size of prey that can be caught
- Throat depth affects how the hook sits in the fish's mouth
- Point angle impacts penetration and retention

### Comparative Analysis
- Compare metrics across different contexts (time periods, regions, etc.)
- Higher variability suggests less standardization in production
- Clusters may represent different functional types or manufacturing traditions

## Advanced Usage

### Customizing Analysis
- Adjust the number of harmonics to control detail level (more harmonics = more detail)
- Increase clusters to identify more subtle shape variations
- Use shape extremes on PCA plot to identify representative specimens

### Batch Export
- The PDF report includes all specimens analyzed in a single session
- Process multiple batches for large collections

## License

This project is open source. Please cite the application if used in academic research.

## Citation

Lipo, Carl (2025). FishhookAnalysis: A Shiny application for morphometric analysis of fishhooks. GitHub repository: https://github.com/clipo/fishhookAnalysis

## Acknowledgments

This application builds upon methods from geometric morphometrics, particularly the Momocs package for R developed by Vincent Bonhomme.