import argparse
import pandas as pd

hrd_cancer_types = ["Breast", "Ovary", "Pancreas", "Prostate", "Fallopian tube"]

# use all hrd samples

# curate training data

# the input sample tsv file can be retrieved from database using this SQL query:

'''
select clinical.sampleId, clinical.hmfPatientId, clinical.primaryTumorLocation, purity.qcStatus, purity.purity, purity.ploidy, chord.BRCA1, chord.BRCA2, chord.hrd, chord.hrStatus, chord.hrdType,
chord.remarksHrStatus, chord.remarksHrdType from clinical,purity,chord
where clinical.sampleId = chord.sampleId and clinical.sampleId = purity.sampleId;
'''
def select_samples(all_samples_tsv, validation_patient_ids, num_extra_random_samples) -> pd.DataFrame:
    df = pd.read_csv(all_samples_tsv, sep="\t")
    df = df[(df["qcStatus"] == "PASS") & (df["hrStatus"] != "CANNOT_BE_DETERMINED")].reset_index(drop=True)

    print(len(df))
    if validation_patient_ids is not None:
        df["validation"] = df["hmfPatientId"].isin(validation_patient_ids)
        if df[df["validation"]].empty:
            print("no sample is in validation set")
            raise RuntimeError("no sample is in validation set")
        df = df[~df["validation"]].drop(columns="validation").reset_index(drop=True)
        print(len(df))

    # we keep all that are the types we are most interested in, and also all HRD deficient ones
    df["mustKeep"] = df["primaryTumorLocation"].isin(hrd_cancer_types) | (df["hrStatus"] == "HR_DEFICIENT")
    print(len(df))
    # the ones we are not designated as keep, we randomly select 500
    df = pd.concat([df[df['mustKeep']], df[~df['mustKeep']].sample(n=num_extra_random_samples, random_state=42)]).drop(columns="mustKeep")
    return df


def main():
    parser = argparse.ArgumentParser(description="train hrd predictor")
    parser.add_argument('--sample_tsv', help='sample tsv file', required=True)
    parser.add_argument('--validation_sample_tsv', help='validation samples tsv file, these samples will be excluded', required=False)
    parser.add_argument('--num_extra_random_samples', help='number of randomly selected extra samples', required=True)
    parser.add_argument('--output_tsv', help='output sample TSV file', required=True)
    args = parser.parse_args()

    validation_samples = None

    if args.validation_sample_tsv:
        validation_samples = pd.read_csv(args.validation_sample_tsv, sep="\t")["hmfPatientId"]

    df = select_samples(args.sample_tsv, validation_samples, args.num_extra_random_samples)
    df.to_csv(args.output_tsv, sep="\t", index=False)

if __name__ == "__main__":
    main()
