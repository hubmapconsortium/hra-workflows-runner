## Anatomogram (EBI Single Cell Expression Atlas)

Datasets are extracted from the `{CODE}.project.h5ad` files on the
[SCEA FTP site](https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments/).
Each experiment is split into one dataset per donor (and sampling site when available).
Dataset ids have the form `ANATOMOGRAM-{CODE}-{individual}[-{sampling_site}]`.

Author assigned cell types (`authors_cell_type_-_ontology_labels`) are summarized using the `author` algorithm.
The corresponding crosswalk rows (`Organ_Level` = `anatomogram`) in `crosswalking-tables/author.csv`
can be regenerated with `python3 src/anatomogram/build_crosswalk.py`.

## Environment configuration options

#### Required
None

#### Optional
- `ANATOMOGRAM_EXPERIMENTS` - Comma separated list of experiment codes to process (default: `E-CURD-119,E-MTAB-10553,E-GEOD-130148,E-MTAB-5061`)
- `ANATOMOGRAM_COUNTS_LAYER` - Layer in the h5ad files containing counts (default: `filtered`)
