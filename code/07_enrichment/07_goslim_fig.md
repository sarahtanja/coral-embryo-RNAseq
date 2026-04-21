# Visualize GOslims by numbers of DEGs

2026-02-06

- [<span class="toc-section-number">1</span> Background](#background)
- [<span class="toc-section-number">2</span> Load inputs](#load-inputs)
- [<span class="toc-section-number">3</span> Wrangle](#wrangle)
- [<span class="toc-section-number">4</span> What is
  unique?](#what-is-unique)
- [<span class="toc-section-number">5</span> What is
  shared?](#what-is-shared)
- [<span class="toc-section-number">6</span> How do I show this in
  context of background universe and stage
  signals?](#how-do-i-show-this-in-context-of-background-universe-and-stage-signals)

# Background

# Load inputs

``` r
input_path <- "../../output/07_enrichment/goslim/"
```

``` r
goslim_summary <- read_csv(file.path(input_path, "goslim_summary.csv"))
```

``` r
glimpse(goslim_summary)
```

``` r
summary(goslim_summary)
```

``` r
colSums(is.na(goslim_summary))
```

# Wrangle

remove rows that are 0 or NA

``` r
goslim_summary %>%
  drop_na(genes) %>% 
  filter(DEdueto == "leachate") %>% 
  filter(n_genes > 2) %>%
  group_by(module) %>%
  slice_max(order_by = n_genes, n = 10) %>%
  ungroup() %>% 
ggplot(., aes(x = n_genes, y = Term, fill = Term)) +
  geom_col(position = position_dodge(width = 0.9), alpha = 0.85) +
  facet_wrap(~ pattern, scales = "free_y", ncol = 1) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(
    title = "GOslim BP Terms: Number of Genes Mapped",
    x = "Number of differentially expressed genes mapped to GOslim terms",
    y = "Top 10 GOslim terms per gene expression module"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank())
```

``` r
ggsave(
  filename = "goslim_summary_barplot.png",
  plot = last_plot(),
  path = "../../output/07_enrichment/goslim/",
  width = 10,
  height = 8,
  units = "in",
  dpi = 900
)
```

# What is unique?

# What is shared?

# How do I show this in context of background universe and stage signals?

How many GOSlims are there?

``` r
length(unique(goslim_summary$Term))
```

Which are the top terms of *leachate-responsive DEGs* represented in
each pattern?

``` r
goslim_summary %>%
  filter(DEdueto == "leachate", n_genes > 2) %>%
  group_by(module) %>%
  slice_max(order_by = Percent, n = 15) %>%
  ungroup() %>%
  count(module, Term) %>%
  count(module, name = "n_unique_terms")
```

``` r
top_terms <- goslim_summary %>%
  filter(DEdueto == "leachate", n_genes > 2) %>%
  group_by(module) %>%
  slice_max(order_by = Percent, n = 15, with_ties = FALSE) %>%
  ungroup()

overlap_terms <- top_terms %>%
  group_by(Term) %>%
  summarise(
    n_modules = n_distinct(module),
    modules = paste(sort(unique(module)), collapse = ", "),
    patterns = paste(sort(unique(pattern)), collapse = ", "),
    mean_percent = mean(Percent, na.rm = TRUE),
    mean_n_genes = mean(n_genes, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_modules > 1) %>%
  arrange(desc(n_modules), desc(mean_percent))

overlap_terms
```

``` r
overlap_detail <- top_terms %>%
  semi_join(overlap_terms, by = "Term") %>%
  arrange(Term, module, desc(Percent)) %>%
  select(Term, module, pattern, Percent, n_genes, Count, GO.IDs)

overlap_detail
```

``` r
goslim_summary %>%
  drop_na(genes) %>% 
  filter(DEdueto == "leachate") %>% 
  filter(n_genes > 2) %>%
  group_by(module) %>%
  slice_max(order_by = Percent, n = 10) %>%
  summarise(
    unique_terms = list(unique(Term)),
    n_unique_terms = n_distinct(Term),
    .groups = "drop"
  )
```

``` r
library(tidyverse)

gos_plot <- goslim_summary %>%
  drop_na(genes) %>%
  filter(DEdueto == "leachate", n_genes > 2)

term_levels <- gos_plot %>%
  group_by(Term) %>%
  summarise(total_n_genes = sum(n_genes, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_n_genes), Term) %>%
  slice_head(n = 10) %>%
  pull(Term)

circ_df <- gos_plot %>%
  filter(Term %in% term_levels) %>%
  group_by(pattern, Term) %>%
  summarise(n_genes = sum(n_genes), .groups = "drop") %>%
  complete(pattern, Term = term_levels, fill = list(n_genes = 0)) %>%
  mutate(
    pattern = factor(pattern, levels = c("down", "hump", "up")),
    Term = factor(Term, levels = term_levels),
    Term_wrap = str_wrap(as.character(Term), width = 18)
  )

ggplot(circ_df, aes(x = factor(Term_wrap, levels = unique(Term_wrap)), y = n_genes, fill = Term)) +
  geom_col(width = 0.95, alpha = 0.9) +
  coord_polar(start = -pi / 2) +
  facet_wrap(~pattern, nrow = 1) +
  labs(
    title = "GOslim BP terms across expression patterns",
    x = NULL,
    y = "Number of DE genes mapped to term"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 8),
    legend.position = "none"
  )
```

``` r
library(tidyverse)

gos_plot <- goslim_summary %>%
  drop_na(genes) %>%
  filter(DEdueto == "leachate", n_genes > 2)

# choose one shared set of terms across all patterns
term_levels <- gos_plot %>%
  group_by(Term) %>%
  summarise(total_percent = sum(Percent, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_percent), Term) %>%
  slice_head(n = 10) %>%
  pull(Term)

circ_df <- gos_plot %>%
  filter(Term %in% term_levels) %>%
  group_by(pattern, Term) %>%
  summarise(
    Percent = max(Percent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  complete(pattern, Term = term_levels, fill = list(Percent = 0)) %>%
  mutate(
    pattern = factor(pattern, levels = c("down", "hump", "up")),
    Term = factor(Term, levels = term_levels),
    Term_wrap = str_wrap(as.character(Term), width = 18)
  )

ggplot(
  circ_df,
  aes(
    x = factor(Term_wrap, levels = unique(Term_wrap)),
    y = Percent,
    fill = Term
  )
) +
  geom_col(width = 0.95, alpha = 0.9) +
  coord_polar(start = -pi / 2) +
  facet_wrap(~pattern, nrow = 1) +
  labs(
    title = "GOslim BP terms across expression patterns",
    subtitle = "Bar length reflects percent representation",
    x = NULL,
    y = "Percent"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 8),
    legend.position = "none"
  )

ggplot(
  circ_df,
  aes(
    x = factor(Term_wrap, levels = unique(Term_wrap)),
    y = Percent,
    fill = Term
  )
) +
  geom_col(width = 0.95, alpha = 0.9) +
  coord_polar(start = -pi / 2) +
  facet_wrap(~pattern, nrow = 1) +
  scale_y_continuous(labels = scales::label_percent(scale = 1)) +
  labs(
    title = "GOslim BP terms across expression patterns",
    subtitle = "Bar length reflects percent representation",
    x = NULL,
    y = "Percent of genes mapped to GOslim term"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 8),
    legend.position = "none"
  )
```
