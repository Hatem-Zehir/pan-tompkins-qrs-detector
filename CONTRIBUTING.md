# Contributing to pantompkin

Thank you for your interest in contributing! Bug reports, feature suggestions,
documentation improvements, and code contributions are all welcome.

---

## Getting Started

### 1. Fork and clone the repository

```bash
git clone https://github.com/Hatem-Zehir/pan-tompkins-qrs-detector.git
cd pan-tompkins-qrs-detector
```

### 2. Set up a development environment

```bash
python -m venv .venv
source .venv/bin/activate        # Linux / macOS
.venv\Scripts\activate           # Windows

pip install -e ".[dev]"
```

### 3. Run the tests

```bash
pytest
```

All tests must pass before you submit a pull request.

---

## How to Contribute

### Reporting bugs

Open an issue on [GitHub](https://github.com/Hatem-Zehir/pan-tompkins-qrs-detector/issues) and include:
- A clear description of the bug
- A minimal reproducible example
- Your Python version and operating system

### Suggesting features

Open a GitHub issue with the label `enhancement`. Describe the feature and
explain why it would be useful to the research community.

### Submitting a pull request

1. Create a branch from `main`:
   ```bash
   git checkout -b fix/your-fix-name
   ```
2. Make your changes and add or update tests as needed.
3. Confirm all tests pass: `pytest`
4. Commit with a clear message and push to your fork.
5. Open a pull request against the `main` branch.

---

## Code Style

- Follow [PEP 8](https://peps.python.org/pep-0008/).
- Write docstrings in NumPy format for all public functions.
- Keep functions focused and independently testable.

---

## Attribution

Contributors will be acknowledged in the `CONTRIBUTORS` file and in release notes.

---

## Code of Conduct

All contributors are expected to follow our [Code of Conduct](CODE_OF_CONDUCT.md).
