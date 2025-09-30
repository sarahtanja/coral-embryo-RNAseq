## Recommendations for Annotation/Enrichment of a Single Gene

### 1. GO enrichment on a single gene is not appropriate
- GO enrichment tests for overrepresentation of GO terms in a set of genes versus a background. With only one gene, there is nothing to compare and no meaningful statistic can be computed.
- Most tools will either error or trivially return the GO terms of your one gene.

**Best practice:**
- If you have only one gene in a cluster, report its functional annotation directly.

### 2. How to summarize and report annotation for a single gene
- Gather all available functional annotation fields for the gene (Gene ID, Description, GO terms, KEGG pathways, COG, PFAMs, etc.).
- Present these in a readable table or narrative. Example:

| Field       | Value                                  |
|-------------|----------------------------------------|
| Gene ID     | Montipora_capitata_HIv3___TS.g15697.t1 |
| Description | regulation of mRNA cap binding         |
| GO Terms    | (list GO terms & names below)          |
| KEGG        | ko03013, ko05416, map03013, map05416   |
| COG         | J                                      |
| PFAMs       | MA3, MIF4G, W2                         |

- **Narrative example:**
  > The gene Montipora_capitata_HIv3___TS.g15697.t1 is described as "regulation of mRNA cap binding". It is annotated with GO terms including mRNA cap binding, nuclear division, reproduction, etc., and is associated with KEGG pathways ko03013, ko05416, etc.

### 3. Map GO term IDs to names in R
```r
# Extract GO terms as character vector
go_ids <- unlist(strsplit(sig_genes_af4$GOs[1], ","))
# Use GO.db to map to names:
library(GO.db)
go_names <- Term(go_ids)
go_map <- data.frame(GO_ID = go_ids, GO_Name = go_names, stringsAsFactors = FALSE)
# View or print:
print(go_map)
```
- You may want to add the Ontology column via `Ontology(go_ids)`.

### 4. Reporting
- For clusters with just one gene, do not attempt GO enrichment—just report the gene’s annotation as described above.

---

Copilot is powered by AI, so mistakes are possible. Leave a comment via the 👍 👎 to share your feedback and help improve the experience.
