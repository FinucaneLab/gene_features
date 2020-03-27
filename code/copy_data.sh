#!/bin/bash

# Copy data from google bucket
gsutil -m rsync -r gs://gene_features/data/ ../data/



