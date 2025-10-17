package com.hartwig.hmftools.lilac.qc;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.FAIL_LOW_COVERAGE_THRESHOLD;
import static com.hartwig.hmftools.lilac.LilacConstants.WARN_INDEL_THRESHOLD;
import static com.hartwig.hmftools.lilac.LilacConstants.WARN_LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.lilac.LilacConstants.WARN_LOW_COVERAGE_THRESHOLD;
import static com.hartwig.hmftools.lilac.LilacConstants.WARN_UNMATCHED_HAPLOTYPE_SUPPORT;
import static com.hartwig.hmftools.lilac.qc.LilacQCStatus.FAIL;
import static com.hartwig.hmftools.lilac.qc.LilacQCStatus.WARN_LOW_BASE_QUAL;
import static com.hartwig.hmftools.lilac.qc.LilacQCStatus.WARN_LOW_COVERAGE;
import static com.hartwig.hmftools.lilac.qc.LilacQCStatus.WARN_UNMATCHED_ALLELE;
import static com.hartwig.hmftools.lilac.qc.LilacQCStatus.WARN_UNMATCHED_AMINO_ACID;
import static com.hartwig.hmftools.lilac.qc.LilacQCStatus.WARN_UNMATCHED_HAPLOTYPE;
import static com.hartwig.hmftools.lilac.qc.LilacQCStatus.WARN_UNMATCHED_INDEL;
import static com.hartwig.hmftools.lilac.qc.LilacQCStatus.WARN_UNMATCHED_SOMATIC_VARIANT;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.hla.LilacQcData;
import com.hartwig.hmftools.common.utils.file.FileLock;
import com.hartwig.hmftools.lilac.GeneSelector;
import com.hartwig.hmftools.lilac.hla.HlaAllele;

public final class LilacQC
{
    private final GeneSelector mGenes;

    public final Set<LilacQCStatus> Status;
    public final double ScoreMargin;
    public final String NextSolutionInfo;
    public final HlaAllele HlaYAllele;
    public final int MedianBaseQual;

    public final AminoAcidQC AminoAcidQC;
    public final BamQC BamQC;
    public final CoverageQC CoverageQC;
    public final HaplotypeQC HaplotypeQC;
    public final SomaticVariantQC SomaticVariantQC;

    public LilacQC(
            final GeneSelector genes, double scoreMargin, final String nextSolutionInfo, int medianBaseQual, final HlaAllele hlaYAllele,
            final AminoAcidQC aminoAcidQC, final BamQC bamQC,
            final CoverageQC coverageQC, final HaplotypeQC haplotypeQC, final SomaticVariantQC somaticVariantQC)
    {
        mGenes = genes;

        ScoreMargin = scoreMargin;
        NextSolutionInfo = nextSolutionInfo;
        MedianBaseQual = medianBaseQual;

        HlaYAllele = hlaYAllele;
        AminoAcidQC = aminoAcidQC;
        BamQC = bamQC;
        CoverageQC = coverageQC;
        HaplotypeQC = haplotypeQC;
        SomaticVariantQC = somaticVariantQC;

        Status = Sets.newHashSet();
        populateStatus();
    }

    public List<String> getHeaderItems()
    {
        List<String> columns = Lists.newArrayList();
        columns.add("Genes");
        columns.add(LilacQcData.FLD_QC_STATUS);
        columns.add("ScoreMargin");
        columns.add("NextSolutionAlleles");
        columns.add("MedianBaseQuality");
        columns.add("HlaYAllele");
        columns.addAll(BamQC.header());
        columns.addAll(CoverageQC.header());
        columns.addAll(AminoAcidQC.header());
        columns.addAll(HaplotypeQC.header());
        columns.addAll(SomaticVariantQC.header());
        return columns;
    }

    public List<String> getBodyItems()
    {
        List<String> columns = Lists.newArrayList();
        columns.add(mGenes.toString());

        StringJoiner sj = new StringJoiner(ITEM_DELIM);
        Status.forEach(x -> sj.add(x.toString()));
        columns.add(sj.toString());

        columns.add(String.format("%.3f", ScoreMargin));
        columns.add(NextSolutionInfo);
        columns.add(String.valueOf(MedianBaseQual));
        columns.add(HlaYAllele != null ? HlaYAllele.toString() : "NONE");
        columns.addAll(BamQC.body());
        columns.addAll(CoverageQC.body());
        columns.addAll(AminoAcidQC.body());
        columns.addAll(HaplotypeQC.body());
        columns.addAll(SomaticVariantQC.body());
        return columns;
    }

    public String header()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        getHeaderItems().forEach(sj::add);
        return sj.toString();
    }

    public String body()
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        getBodyItems().forEach(sj::add);
        return sj.toString();
    }

    public void log(final String sampleId)
    {
        // LL_LOGGER.info("QC Stats:");
        StringJoiner sj = new StringJoiner(",");
        List<String> headerItems = getHeaderItems();
        List<String> bodyItems = getBodyItems();

        for(int i = 0; i < headerItems.size(); ++i)
        {
            sj.add(String.format("%s=%s", headerItems.get(i), bodyItems.get(i)));
        }

        LL_LOGGER.info("{} QC Stats: {}", sampleId, sj.toString());
    }

    public void writefile(final String fileName)
    {
        File file = new File(fileName);
        List<String> existingRecords = Lists.newArrayList();
        try(FileLock fileLock = FileLock.create(file))
        {
            BufferedReader reader = fileLock.getBufferedReader();
            reader.readLine();
            String line = reader.readLine();
            while(line != null)
            {
                String genesField = line.split(TSV_DELIM)[0];
                GeneSelector genes;
                try
                {
                    genes = GeneSelector.valueOf(genesField);
                }
                catch(IllegalArgumentException e)
                {
                    genes = null;
                }

                if(genes != null && genes != mGenes)
                    existingRecords.add(line);

                line = reader.readLine();
            }

            fileLock.clear();
            BufferedWriter writer = fileLock.getBufferedWriter();
            writer.write(header());
            writer.newLine();

            for(String existingRecord : existingRecords)
            {
                writer.write(existingRecord);
                writer.newLine();
            }

            writer.write(body());
            writer.newLine();
            writer.flush();
        }
        catch(Exception e)
        {
            LL_LOGGER.error("failed to update {}: {}", fileName, e.toString());
        }
    }

    private void populateStatus()
    {
        if(CoverageQC.TotalFragments == 0 || BamQC.totalLowCoverage() >= FAIL_LOW_COVERAGE_THRESHOLD)
        {
            Status.add(FAIL);
            return;
        }

        if(MedianBaseQual < WARN_LOW_BASE_QUAL_THRESHOLD)
        {
            Status.add(WARN_LOW_BASE_QUAL);
        }

        if(BamQC.totalLowCoverage() >= WARN_LOW_COVERAGE_THRESHOLD)
        {
            Status.add(WARN_LOW_COVERAGE);
        }

        double haplotypeWarnThreshold = CoverageQC.TotalFragments * WARN_UNMATCHED_HAPLOTYPE_SUPPORT;
        if(HaplotypeQC.UnusedHaplotypeMaxFrags >= haplotypeWarnThreshold)
        {
            Status.add(WARN_UNMATCHED_HAPLOTYPE);
        }

        if(AminoAcidQC.UnusedAminoAcidMaxFrags >= haplotypeWarnThreshold)
        {
            Status.add(WARN_UNMATCHED_AMINO_ACID);
        }

        if(CoverageQC.CountsByGene.values().stream().anyMatch(x -> x == 0))
        {
            Status.add(WARN_UNMATCHED_ALLELE);
        }

        double indelWarnThreshold = CoverageQC.TotalFragments * WARN_INDEL_THRESHOLD;
        if(BamQC.DiscardedIndelMaxFrags >= indelWarnThreshold)
        {
            Status.add(WARN_UNMATCHED_INDEL);
        }

        if(SomaticVariantQC.unmatchedVariants() > 0)
        {
            Status.add(WARN_UNMATCHED_SOMATIC_VARIANT);
        }

        if(Status.isEmpty())
        {
            Status.add(LilacQCStatus.PASS);
        }
    }
}
