# Metabolome preprocessing
n_met_feats <- ncol(raw.metabolome) - 3
met_feat_names <- paste0("met.feature", seq_len(n_met_feats))
metabolome <- raw.metabolome
colnames(metabolome)[4:ncol(metabolome)] <- met_feat_names
metabolome[, 4:ncol(metabolome)] <- scale(metabolome[, 4:ncol(metabolome)])
map.metabolome <- data.frame(name = colnames(raw.metabolome)[4:ncol(raw.metabolome)], index = met_feat_names)

# Microbiome preprocessing
n_mb_feats <- ncol(raw.microbiome) - 3
mb_feat_names <- paste0("mb.feature", seq_len(n_mb_feats))
microbiome <- raw.microbiome
colnames(microbiome)[4:ncol(microbiome)] <- mb_feat_names
feature_cols <- mb_feat_names
non_zero_cols <- feature_cols[colSums(microbiome[, feature_cols]) != 0]
microbiome <- microbiome[, c(colnames(microbiome)[1:3], non_zero_cols)]
microbiome[, 4:ncol(microbiome)] <- scale(microbiome[, 4:ncol(microbiome)])
map.microbiome <- data.frame(name = colnames(raw.microbiome)[4:ncol(raw.microbiome)], index = mb_feat_names)
