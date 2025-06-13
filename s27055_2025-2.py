from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt


def search_taxid(email, api_key, taxid):
    Entrez.email = email
    Entrez.api_key = api_key
    try:
        handle = Entrez.esearch(db="nucleotide", term=f"txid{taxid}[Organism]", usehistory="y")
        search_results = Entrez.read(handle)
        count = int(search_results["Count"])
        if count == 0:
            print(f"No records found for TaxID {taxid}.")
            return None, None, None
        return count, search_results["WebEnv"], search_results["QueryKey"]
    except Exception as e:
        print(f"Error searching TaxID {taxid}: {e}")
        return None, None, None


def fetch_records(webenv, query_key, max_records):
    try:
        handle = Entrez.efetch(
            db="nucleotide",
            rettype="gb",
            retmode="text",
            retmax=max_records,
            webenv=webenv,
            query_key=query_key
        )
        return SeqIO.parse(handle, "genbank")
    except Exception as e:
        print(f"Error fetching records: {e}")
        return []


def main():
    email = input("Podaj swój adres e-mail: ")
    api_key = input("Podaj swój klucz API (może być pusty): ")
    taxid = input("Podaj identyfikator taksonomiczny (TaxID): ")
    min_length = int(input("Podaj minimalną długość sekwencji: "))
    max_length = int(input("Podaj maksymalną długość sekwencji: "))
    csv_filename = input("Podaj nazwę pliku CSV (np. wyniki.csv): ")
    plot_filename = input("Podaj nazwę pliku wykresu (np. wykres.png): ")

    count, webenv, query_key = search_taxid(email, api_key, taxid)
    if not count:
        return

    print(f"Found {count} records. Fetching data...")
    records = fetch_records(webenv, query_key, max_records=count)

    filtered_data = []
    for record in records:
        seq_length = len(record.seq)
        if min_length <= seq_length <= max_length:
            filtered_data.append((record.id, seq_length, record.description))

    if not filtered_data:
        print("No sequences found within the specified length range.")
        return

    df = pd.DataFrame(filtered_data, columns=["Accession", "Length", "Description"])
    df.to_csv(csv_filename, index=False)
    print(f"Saved filtered data to {csv_filename}.")

    df = df.sort_values("Length", ascending=False)
    plt.figure(figsize=(10, 6))
    plt.plot(df["Accession"], df["Length"], marker="o")
    plt.xticks(rotation=90)
    plt.xlabel("Accession")
    plt.ylabel("Length")
    plt.title("Sequence Lengths")
    plt.tight_layout()
    plt.savefig(plot_filename)
    print(f"Saved plot to {plot_filename}.")


if __name__ == "__main__":
    main()
