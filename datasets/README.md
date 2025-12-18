# Datasets

Scalc loads datasets from CSV files. Each file can have:
- Comment lines starting with `#` (displayed as metadata)
- A header row with column names
- Data rows (numeric and/or text columns)

## Usage

```
>>> load iris           # Load datasets/iris.csv
>>> load hald           # Load datasets/hald.csv
>>> load mydata.csv     # Load any CSV file
```

Data is loaded into:
- `data`: Numeric columns as a matrix
- `labels`: Text/categorical column encoded as integers (1, 2, 3...)

## Built-in Datasets

### iris.csv - Fisher's Iris Dataset
150 samples, 4 features, 3 classes. Classic classification benchmark.
- **Columns**: sepal_length, sepal_width, petal_length, petal_width, species
- **Classes**: setosa (1), versicolor (2), virginica (3)

### hald.csv - Hald Cement Dataset  
13 samples, 4 predictors, 1 response. Classic regression benchmark.
- **Predictors**: calcium_aluminate, dicalcite_silicate, calcium_aluminoferrite, beta_dicalcium_silicate
- **Response**: heat (calories/gram evolved during hardening)

### imports85.csv - UCI Automobile Dataset
109 samples (complete cases), 16 numeric features.
- Missing values marked with `?` are loaded as NaN
- **Features**: symboling, normalized_losses, wheel_base, length, width, height, curb_weight, engine_size, bore, stroke, compression_ratio, horsepower, peak_rpm, city_mpg, highway_mpg, price

## CSV Format

```csv
# Optional comment line (metadata)
# Another comment
column1,column2,column3
1.0,2.0,category_a
3.0,4.0,category_b
```

- Missing values: Use `?` for missing numeric data (loaded as NaN)
- Text columns: Automatically converted to numeric codes (1, 2, 3...)

## Environment Variable

Set `SC_DATASETS` to specify an alternate dataset directory:
```bash
export SC_DATASETS=/path/to/datasets
sc
>>> load mydata
```
