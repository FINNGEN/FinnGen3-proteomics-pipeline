# Installation Instructions

## System Requirements

- **R**: Version 4.0 or higher
- **Operating System**: Linux (recommended), macOS, or Windows (with WSL)
- **Memory**: Minimum 8GB RAM (16GB+ recommended)
- **Storage**: Sufficient space for input data and outputs

## Installation Methods

### Method 1: Docker (Recommended)

Docker provides a consistent environment with all dependencies pre-installed.

#### Prerequisites
- Docker Engine 20.10+
- Docker Compose 2.0+ (optional, for docker-compose)

#### Steps

1. **Clone the repository**
```bash
git clone <repository-url>
cd fg3_olink_pipeline
```

2. **Build the Docker image**
```bash
docker build -t fg3-olink-pipeline .
```

3. **Verify installation**
```bash
docker run fg3-olink-pipeline Rscript --version
```

### Method 2: Local R Installation

#### Prerequisites
- R 4.0+ installed
- R development tools (build-essential on Linux, Xcode on macOS)

#### Steps

1. **Clone the repository**
```bash
git clone <repository-url>
cd fg3_olink_pipeline
```

2. **Install R packages**
```bash
Rscript install_packages.R
```

3. **Install BioConductor packages**
```bash
R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('sva')"
```

4. **Install OlinkAnalyze** (if needed)
```bash
R -e "install.packages('OlinkAnalyze', repos='https://cloud.r-project.org/')"
```

5. **Install external tools**

   **PLINK** (required for pQTL steps):
   ```bash
   # Linux
   wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip
   unzip plink_linux_x86_64_20230116.zip
   sudo mv plink /usr/local/bin/

   # macOS
   brew install plink
   ```

   **Google Cloud SDK** (required for GCS access in pQTL steps):
   ```bash
   # Follow instructions at: https://cloud.google.com/sdk/docs/install
   ```

6. **Verify installation**
```bash
Rscript -e "library(data.table); library(arrow); library(yaml); cat('Core packages loaded successfully!\n')"
```

### Method 3: Conda Environment

1. **Create conda environment**
```bash
conda create -n fg3-olink-pipeline r-base=4.3 r-essentials
conda activate fg3-olink-pipeline
```

2. **Install R packages**
```bash
Rscript install_packages.R
```

3. **Install external tools** (see Method 2, step 5)

## Configuration

1. **Copy config template**
```bash
cp config/config.yaml.template config/config.yaml
```

2. **Edit config.yaml**
   - Set `data.npx_matrix_file` to your pre-filtered NPX matrix path
   - Set `data.metadata_file` to your metadata file path
   - Configure other paths and parameters as needed

3. **Set environment variables** (optional)
```bash
export PIPELINE_CONFIG=/path/to/config.yaml
export PIPELINE_BATCH_ID=batch_01
export PIPELINE_OUTPUT_DIR=/path/to/output
```

## Testing the Installation

### Test with Docker

```bash
docker run fg3-olink-pipeline Rscript scripts/run_pipeline.R --help
```

### Test Local Installation

```bash
# Test package loading
Rscript -e "
library(data.table)
library(tidyverse)
library(arrow)
library(yaml)
library(logger)
cat('All core packages loaded successfully!\n')
"

# Test script execution
export PIPELINE_CONFIG=config/config.yaml
Rscript scripts/00_data_loader.R --help
```

## Optional Dependencies

### Machine Learning Models (for sex prediction)

The sex outlier detection step can use optional ML models:
- **Random Forest**: `install.packages('randomForest')`
- **XGBoost**: `install.packages('xgboost')`
- **Keras/TensorFlow**: `install.packages('keras3'); keras3::install_tensorflow()`

These are optional - the pipeline will use available models or fall back to elastic-net.

### Google Cloud Access (for pQTL steps)

If using pQTL outlier detection with GCS data:

1. **Install Google Cloud SDK**
2. **Authenticate**
```bash
gcloud auth login
gcloud auth application-default login
```

3. **Set credentials in environment**
```bash
export GOOGLE_APPLICATION_CREDENTIALS=/path/to/credentials.json
```

## Troubleshooting

### R Package Installation Issues

**Problem**: Package installation fails
- **Solution**: Update R to latest version
- **Solution**: Install development tools (build-essential, libcurl-dev, etc.)

**Problem**: BioConductor packages fail
- **Solution**: Update BiocManager: `BiocManager::install(version = "3.18")`

### Docker Issues

**Problem**: Docker build fails
- **Solution**: Ensure Docker has sufficient resources (memory, disk space)
- **Solution**: Check internet connection for package downloads

**Problem**: Permission denied in Docker volumes
- **Solution**: Check file permissions on host directories
- **Solution**: Use `--user` flag to run as specific user

### External Tool Issues

**Problem**: PLINK not found
- **Solution**: Ensure PLINK is in PATH or specify full path in config
- **Solution**: Verify PLINK installation: `plink --version`

**Problem**: GCS access fails
- **Solution**: Verify Google Cloud SDK is installed and authenticated
- **Solution**: Check `GOOGLE_APPLICATION_CREDENTIALS` environment variable

## Next Steps

After installation:

1. **Prepare your data**: Place pre-filtered NPX matrix and metadata in `data/`
2. **Configure pipeline**: Edit `config/config.yaml`
3. **Run pipeline**: See [README.md](README.md) for usage instructions

## Support

For installation issues:
- Check the troubleshooting section above
- Review error messages carefully
- Open an issue on GitHub with:
  - Operating system and version
  - R version (`R --version`)
  - Error messages and logs
  - Steps to reproduce




