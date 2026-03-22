# Documentation Site

## Why MkDocs + Material?

This documentation site is designed to feel closer to a modern conference project page than a plain API manual.

The stack is:

- MkDocs
- Material for MkDocs
- custom CSS for the Ariadne visual identity
- Read the Docs deployment configuration

## Local preview

Install documentation dependencies:

```bash
python -m pip install -r docs/requirements.txt
```

Start a local dev server:

```bash
mkdocs serve
```

The site will be available locally at:

```text
http://127.0.0.1:8000/
```

## Production build

To build the static site locally:

```bash
mkdocs build
```

The generated site will be written to:

```text
site/
```

## Read the Docs configuration

The repository includes:

- `.readthedocs.yaml`
- `mkdocs.yml`
- `docs/requirements.txt`

This is enough for a standard Read the Docs build once the repository is connected to an RTD project.

## Site structure

The English documentation site currently includes:

- `Home`
- `Getting Started`
- `Method`
- `Tutorials`
- `CLI Reference`
- `Outputs`
- `Documentation Site`
- `FAQ`
- `Citation`

## Design goals

This documentation site intentionally emphasizes:

- clean project-page aesthetics
- stage-by-stage method explanation
- figure-friendly layout
- practical command examples
- explicit result interpretation guidance

## Suggested future extensions

If you want to evolve this documentation further, strong next additions would be:

- a gallery page with real `embedding.svg` and tree previews
- auto-generated API documentation for internal Python helpers
- dataset preparation guidelines for building new `tree/` reference collections
- versioned documentation releases
