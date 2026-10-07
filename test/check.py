import argparse
from pathlib import Path
import pandas as pd
import json

parser = argparse.ArgumentParser()
parser.add_argument("outputs")
parser.add_argument("name")
args = parser.parse_args()

outputs = Path(args.outputs)
alien_pred = []

acc_to_name = {}

with open('test/data.json', 'r') as f:
    data = json.load(f)

for name, value in data.items():
    acc_to_name[value['accession']] = name

res = {}

for item in outputs.iterdir():
    if item.is_dir() and item.name[:3] == "GCF":
        if not (item / "APP_Alien_genes.tsv").is_file():
            continue

        alien_genes_df = pd.read_csv(
            item / "APP_Alien_genes.tsv",
            sep='\t',
            header=None,
            names=["gene", "type", "age"]
        )

        alien_genes_df = alien_genes_df[
            ~alien_genes_df["gene"].str.startswith("lcl")
        ]

        genes = alien_genes_df["gene"].apply(
            lambda x: f"{x.split(':')[0]}|{item.name}"
        ).to_list()

        alien_pred = [gene for gene in genes if gene.startswith("WP")]

        print(alien_pred[0])

        alien_true = []
        with open(f"test/18genomedb/{item.stem}.1_alien.faa", 'r') as f:
            for line in f:
                if line.startswith('>WP'):
                    alien_true.append(line[1:].strip())

        N = 0
        with open(f"test/18genomedb/{item.stem}.1.faa", 'r') as f:
            for line in f:
                if line.startswith('>'):
                    N += 1

        print(item.stem)

        not_found = list(set(alien_true) - set(alien_pred))
        print("not found:", len(not_found))

        incorrect = list(set(alien_pred) - set(alien_true))
        print("incorrect:", len(incorrect))

        correct = list(set(alien_true) & set(alien_pred))
        print("correct:", len(correct))

        true_aliens = set(alien_true)
        pred_aliens = set(alien_pred)

        TP = len(true_aliens & pred_aliens)
        FP = len(pred_aliens - true_aliens)
        FN = len(true_aliens - pred_aliens)
        TN = N - TP - FP - FN

        accuracy = (TP + TN) / N
        precision = TP / (TP + FP)
        recall = TP / (TP + FN)
        specificity = TN / (TN + FP)
        f1 = 2 * TP / (2 * TP + FP + FN)

        print(
            f"Accuracy: {accuracy}\n"
            f"Precision: {precision}\n"
            f"Recall: {recall}\n"
            f"Specificity: {specificity}\n"
            f"F1: {f1}"
        )

        res[acc_to_name[f"{item.stem}.1"]] = {
            'accession': item.name,
            'accuracy': accuracy,
            'recall': recall,
            'precision': precision,
            'f1': f1,
            'specificity': specificity
        }

        print()

# Create dataframe
df = pd.DataFrame.from_dict(res, orient='index')

# Save individual results
df.to_csv(f'{args.name}.csv')

# Calculate averages
statistics = [
    'accuracy',
    'precision',
    'recall',
    'specificity',
    'f1'
]

averages = df[statistics].mean()

print("=" * 40)
print("AVERAGES")
print("=" * 40)

for statistic in statistics:
    print(f"{statistic.capitalize()}: {averages[statistic]:.6f}")

print("=" * 40)
