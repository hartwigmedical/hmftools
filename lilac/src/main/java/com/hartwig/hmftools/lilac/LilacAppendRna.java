package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.APP_NAME;
import static com.hartwig.hmftools.lilac.ReferenceData.GENE_CACHE;
import static com.hartwig.hmftools.lilac.ReferenceData.NUC_GENE_FRAG_ENRICHMENT;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.lilac.coverage.ComplexBuilder;
import com.hartwig.hmftools.lilac.coverage.ComplexCoverage;
import com.hartwig.hmftools.lilac.coverage.FragmentAlleleMapper;
import com.hartwig.hmftools.lilac.coverage.FragmentAlleles;
import com.hartwig.hmftools.lilac.fragment.AminoAcidFragmentPipeline;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.fragment.NucleotideFragmentFactory;
import com.hartwig.hmftools.lilac.fragment.NucleotideGeneEnrichment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaGene;
import com.hartwig.hmftools.lilac.read.BamRecordReader;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;

import org.jetbrains.annotations.NotNull;

public class LilacAppendRna
{
    private final LilacConfig mConfig;
    private final ReferenceData mRefData;

    public LilacAppendRna(final LilacConfig config, final ConfigBuilder configBuilder)
    {
        mConfig = config;
        mRefData = new ReferenceData(mConfig.ResourceDir, mConfig);
    }

    private static boolean matchesSolutionAllele(final List<LilacAllele> solutionAlleles, final String allele)
    {
        return solutionAlleles.stream().anyMatch(x -> x.allele().equals(allele));
    }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();

        if(!mRefData.load())
        {
            LL_LOGGER.error("reference data loading failed");
            System.exit(1);
        }

        String filename = LilacAllele.generateFilename(mConfig.OutputDir, mConfig.Sample);
        List<LilacAllele> solutionAlleles = Lists.newArrayList();

        try
        {
            solutionAlleles.addAll(LilacAllele.read(filename));
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to load Lilac solution file: {}", filename, e.toString());
            System.exit(1);
        }

        List<HlaSequenceLoci> aminoAcidSequences = Lists.newArrayList();
        List<HlaSequenceLoci> nucleotideSequences = Lists.newArrayList();

        for(LilacAllele allele : solutionAlleles)
        {
            List<HlaSequenceLoci> alleleAaSequences = mRefData.AminoAcidSequences.stream()
                    .filter(x -> allele.allele().equals(x.Allele.toString())).collect(Collectors.toList());

            List<HlaSequenceLoci> alleleNucSequences = mRefData.NucleotideSequences.stream()
                    .filter(x -> allele.allele().equals(x.Allele.asFourDigit().toString())).collect(Collectors.toList());

            if(alleleNucSequences.isEmpty() || alleleAaSequences.isEmpty())
            {
                LL_LOGGER.error("failed to match solution allele({}) to reference data", allele.allele());
                System.exit(1);
            }

            aminoAcidSequences.addAll(alleleAaSequences);
            nucleotideSequences.addAll(alleleNucSequences);
        }

        List<HlaAllele> winningAlleles = aminoAcidSequences.stream().map(x -> x.Allele).collect(Collectors.toList());

        NucleotideFragmentFactory nucleotideFragFactory = new NucleotideFragmentFactory(mRefData);

        AminoAcidFragmentPipeline aminoAcidPipeline = new AminoAcidFragmentPipeline(Collections.emptyList());

        int totalFragmentCount = solutionAlleles.stream().mapToInt(x -> x.refFragments()).sum();

        // recoveredSequences
        // Map<String,Map<Integer,Set<String>>> geneAminoAcidHetLociMap =
                // extractHeterozygousLociSequences(aminoAcidPipeline.getReferenceAminoAcidCounts(), minEvidence, Collections.emptyList());

        Map<HlaGene, Map<Integer, Set<String>>> geneAminoAcidHetLociMap = Maps.newHashMap();

        // Map<String,List<Integer>> refNucleotideHetLociMap = calcNucleotideHeterogygousLoci(mRefNucleotideCounts.heterozygousLoci());
        Map<HlaGene, List<Integer>> refNucleotideHetLociMap = Maps.newHashMap();

        FragmentAlleleMapper fragAlleleMapper = new FragmentAlleleMapper(
                geneAminoAcidHetLociMap, refNucleotideHetLociMap, aminoAcidPipeline.getReferenceNucleotides());

        ComplexCoverage rnaCoverage = extractRnaCoverage(
                mConfig.RnaBam, mConfig, mRefData, nucleotideFragFactory, NUC_GENE_FRAG_ENRICHMENT, aminoAcidPipeline,
                fragAlleleMapper, winningAlleles, aminoAcidSequences, nucleotideSequences);

        /*
        ComplexCoverage refCoverage = !mConfig.tumorOnly() ?
                winningRefCoverage : ComplexCoverage.create(Lists.newArrayList());

        solutionSummary.write(LilacAllele.generateFilename(mConfig.OutputDir, mConfig.Sample));

        mSolutionSummary = SolutionSummary.create(refCoverage, mTumorCoverage, mTumorCopyNumber, mSomaticCodingCounts, mRnaCoverage);
        */

        LL_LOGGER.info("Lilac RNA append complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    public static ComplexCoverage extractRnaCoverage(
            final String rnaBam, final LilacConfig config, final ReferenceData referenceData,
            final NucleotideFragmentFactory nucleotideFragFactory, final NucleotideGeneEnrichment nucleotideGeneEnrichment,
            final AminoAcidFragmentPipeline aminoAcidPipeline, final FragmentAlleleMapper fragAlleleMapper,
            final List<HlaAllele> winningAlleles, final List<HlaSequenceLoci> winningSequences,
            final List<HlaSequenceLoci> winningNucSequences)
    {
        if(rnaBam.isEmpty())
        {
            return ComplexCoverage.create(Lists.newArrayList());
        }

        BamRecordReader rnaBamReader = new BamRecordReader(rnaBam, config, GENE_CACHE.GeneTranscriptMap, nucleotideFragFactory);

        List<Fragment> rnaNucleotideFrags = rnaBamReader.findGeneFragments();

        if(nucleotideGeneEnrichment != null)
            nucleotideGeneEnrichment.checkAddAdditionalGenes(rnaNucleotideFrags);

        List<Fragment> rnaFragments = aminoAcidPipeline.calcComparisonCoverageFragments(rnaNucleotideFrags);

        // mResultsWriter.writeFragments(RNA, rnaFragments);

        LL_LOGGER.info("calculating RNA coverage from frags({} highQual={})", rnaNucleotideFrags.size(), rnaFragments.size());

        List<FragmentAlleles> rnaFragAlleles = fragAlleleMapper.createFragmentAlleles(rnaFragments, winningSequences, winningNucSequences);

        // if(mHlaYCoverage.exceedsThreshold())
        //    mHlaYCoverage.assignFragments(winningAlleles, rnaFragAlleles, rnaFragments, RNA);

        ComplexCoverage rnaCoverage = ComplexBuilder.calcProteinCoverage(rnaFragAlleles, winningAlleles);

        rnaCoverage.populateMissingCoverage(Sets.newLinkedHashSet(winningAlleles));
        rnaCoverage.expandToSixAlleles();

        return rnaCoverage;
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        LilacConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        LilacAppendRna lilacAppendRna = new LilacAppendRna(new LilacConfig(configBuilder), configBuilder);
        lilacAppendRna.run();
    }
}
