from fadns_model.variants import list_variants

for slug, name in list_variants().items():
    print(f"{slug}: {name}")
