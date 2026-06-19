from pathlib import Path
from fadns_model.data_io import load_fatty_acid_csv

for path in sorted(Path("data").glob("*.csv")):
    df = load_fatty_acid_csv(path)
    print(path.name, df.shape, list(df.columns))
