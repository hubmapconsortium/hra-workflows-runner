#!/bin/bash

HRA_GRAPH=https://cdn.humanatlas.io/digital-objects/graph/ccf/latest/graph.ttl
OUTPUT=crosswalks-report.csv

# Convert crosswalks csv to JSON-LD
node src/crosswalks-to-jsonld.js

# Run query to report on number of cell types per organ per annotation tool
comunica-sparql-file -t text/csv -f src/crosswalks-report.rq ./crosswalking-tables/crosswalks.jsonld $HRA_GRAPH > $OUTPUT
