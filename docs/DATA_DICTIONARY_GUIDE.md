# Reading the data dictionary

metadata/file_catalog.csv describes each delimited source table and supplementary table. metadata/data_dictionary.csv contains one row per column, indexed by filename and one-based column position, so blank or repeated source headers remain distinguishable. Observed types, missing-token counts and a first nonmissing example are inferred from a full scan; they are descriptive, not a schema imposed on the data. Numerical values are not changed.

Definitions and units are supplied where they are established by the variable name and assay context. Other fields explicitly retain a source-defined meaning and refer to the corresponding Methods or source header. A generic effect, value, mean or score must not be treated as a common measurement scale across assays. Relative percentages, percentage-point contrasts, log-transformed expression, compensated fluorescence, image-coordinate estimates and program scores have different units.

Original empty strings and tokens such as NA, NaN, nan, NULL and None are counted as missing; zeros remain observed values. Boolean and categorical labels are retained. Matrix column names containing sample identifiers designate columns of a matrix, not independent endpoints. Duplicate or blank headers, if present, are reported in file_catalog.csv and are preserved.

Current canonical tables take precedence over older supporting component exports for interpreting manuscript panel labels. docs/FIGURE_DATA_MAP.csv identifies those current locations and where only partial numerical coverage is provided.
