import pandas as pd
from chemblmltools import chembl_activity_target

df1 = chembl_activity_target(
        db_user='arnau',
        db_password='2405',
        db_name='chembl_35',
        organism_contains='enterobacter',
        max_heavy_atoms=100)

print(df1.head(5))