# Storm

A high-performance data processing library with both Rust and Python interfaces.

## Requirements

- Rust 1.87.0 or later
- Python 3.7 or later
- Cargo (Rust package manager)
- pip (Python package manager)
- Docker (optional, for containerized builds)

## Installation

### Rust Binary

To install the Rust binary:

```bash
# Clone the repository
git clone https://github.com/broadinstitute/storm.git
cd storm

# Build the release binary
cargo build --release --features python

# The binary will be available at target/release/storm
```

### Python Package

To install the Python package:

```bash
# Create and activate a virtual environment
python -m venv venv
source venv/bin/activate  # On Windows, use: venv\Scripts\activate

# Install development dependencies
pip install -r dev-requirements.txt

# Build and install the package
maturin build --release
pip install target/wheels/*.whl
```

### Docker

To build and run using Docker:

```bash
# Build the Docker image
docker build -t storm ./docker

# Run the container
docker run storm
```

## Development

### Setting up the development environment

1. Install Rust toolchain:
   ```bash
   rustup install stable
   rustup default stable
   ```

2. Set up Python environment:
   ```bash
   # Create and activate a virtual environment
   python -m venv venv
   source venv/bin/activate  # On Windows, use: venv\Scripts\activate

   # Install development dependencies
   pip install -r dev-requirements.txt
   ```

### Building from source

1. Build the Rust binary:
   ```bash
   cargo build --release --features python
   ```

2. Build the Python package:
   ```bash
   maturin build --release
   ```

### Running tests

```bash
# Run Rust tests
cargo test --verbose --features python

# Run Python tests
cd python
pytest
```

### Code formatting

```bash
# Format Rust code
cargo fmt

# Format Python code
cd python
black .
isort .
```

## License

[Add your license information here]

## Contributing

[Add contribution guidelines here]