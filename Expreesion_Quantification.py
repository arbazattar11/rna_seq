import pandas as pd

def quantify_expression(gtf_file):
    # Load the GTF file into a DataFrame
    columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    df = pd.read_csv(gtf_file, sep='\t', header=None, comment='#', names=columns)

    # Filter rows to retain only transcripts (feature == 'transcript')
    transcripts = df[df['feature'] == 'transcript']

    # Extract gene_id and transcript_id from attributes column
    transcripts['gene_id'] = transcripts['attributes'].str.extract(r'gene_id "(.*?)";')
    transcripts['transcript_id'] = transcripts['attributes'].str.extract(r'transcript_id "(.*?)";')

    # Group by gene_id and sum up the length of transcripts to get expression
    expression = transcripts.groupby('gene_id')['end'].sum()

    return expression

if __name__ == "__main__":
    # Path to the output GTF file from StringTie
    gtf_file = 'path/to/your/output.gtf'

    # Perform expression quantification
    expression = quantify_expression(gtf_file)

    # Save expression values to a CSV file
    expression.to_csv('expression.csv')
