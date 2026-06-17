# Documentation Site

## Stack

This documentation is built to read like a modern conference project page rather than a plain API manual:

- [MkDocs](https://www.mkdocs.org/)
- [Material for MkDocs](https://squidfunk.github.io/mkdocs-material/)
- custom CSS for the Ariadne visual identity (`docs/stylesheets/extra.css`)
- Read the Docs deployment configuration

## Local preview

```bash
python -m pip install -r docs/requirements.txt
mkdocs serve     # serves at http://127.0.0.1:8000/
```

## Production build

```bash
mkdocs build     # static site written to site/
```

## Read the Docs

The repository includes everything Read the Docs needs once the project is connected:

- `.readthedocs.yaml`
- `mkdocs.yml`
- `docs/requirements.txt`

## Site structure

| Page | Purpose |
|---|---|
| Home | overview and quick start |
| Getting Started | installation and first run |
| Method | the four-stage design, stage by stage |
| ESM Type Model | the ESM2 CeeSs classifier |
| Tutorials | practical analysis pathways |
| CLI Reference | every command and parameter |
| Outputs | every artifact the pipeline produces |
| Documentation Site | this page |
| FAQ | common questions |
| Citation | how to cite Ariadne |

## Design goals

The site deliberately emphasizes clean project-page aesthetics, a stage-by-stage method narrative, a figure-friendly layout, concrete command examples, and explicit guidance on interpreting results.

## Possible extensions

- a gallery page with real `embedding.svg` and tree previews;
- auto-generated API documentation for the Python helpers;
- guidelines for building new `tree/` reference collections;
- versioned documentation releases.
