package com.hartwig.hmftools.lilac;

import static java.lang.Math.floor;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.APP_NAME;
import static com.hartwig.hmftools.lilac.LilacConstants.BASE_QUAL_PERCENTILE;
import static com.hartwig.hmftools.lilac.LilacConstants.LOW_BASE_QUAL_THRESHOLD;
import static com.hartwig.hmftools.lilac.LilacConstants.MIN_EVIDENCE_FACTOR;
import static com.hartwig.hmftools.lilac.LilacConstants.WARN_LOW_COVERAGE_DEPTH;
import static com.hartwig.hmftools.lilac.ReferenceData.GENE_CACHE;
import static com.hartwig.hmftools.lilac.ReferenceData.HLA_CONTEXT_FACTORY;
import static com.hartwig.hmftools.lilac.ReferenceData.NUC_GENE_FRAG_ENRICHMENT;
import static com.hartwig.hmftools.lilac.coverage.HlaComplex.findDuplicates;
import static com.hartwig.hmftools.lilac.evidence.NucleotideFiltering.calcNucleotideHeterogygousLoci;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.CANDIDATE;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.SOLUTION;
import static com.hartwig.hmftools.lilac.fragment.FragmentSource.TUMOR;
import static com.hartwig.hmftools.lilac.fragment.NucleotideFragmentFactory.calculateGeneCoverage;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_A;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_B;
import static com.hartwig.hmftools.lilac.hla.HlaGene.HLA_C;
import static com.hartwig.hmftools.lilac.seq.SequenceCount.extractHeterozygousLociSequences;
import static com.hartwig.hmftools.lilac.variant.SomaticCodingCount.addVariant;

import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.perf.StackSampler;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.lilac.coverage.AlleleCoverage;
import com.hartwig.hmftools.lilac.coverage.ComplexBuilder;
import com.hartwig.hmftools.lilac.coverage.ComplexCoverage;
import com.hartwig.hmftools.lilac.coverage.ComplexCoverageCalculator;
import com.hartwig.hmftools.lilac.coverage.ComplexCoverageRanking;
import com.hartwig.hmftools.lilac.coverage.FragmentAlleleMapper;
import com.hartwig.hmftools.lilac.coverage.FragmentAlleles;
import com.hartwig.hmftools.lilac.coverage.HlaComplex;
import com.hartwig.hmftools.lilac.coverage.HlaComplexFile;
import com.hartwig.hmftools.lilac.coverage.HlaYCoverage;
import com.hartwig.hmftools.lilac.evidence.Candidates;
import com.hartwig.hmftools.lilac.evidence.PhasedEvidence;
import com.hartwig.hmftools.lilac.fragment.AminoAcidFragmentPipeline;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.fragment.NucleotideFragmentFactory;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaGene;
import com.hartwig.hmftools.lilac.qc.AminoAcidQC;
import com.hartwig.hmftools.lilac.qc.BamQC;
import com.hartwig.hmftools.lilac.qc.CoverageQC;
import com.hartwig.hmftools.lilac.qc.HaplotypeQC;
import com.hartwig.hmftools.lilac.qc.LilacQC;
import com.hartwig.hmftools.lilac.qc.SolutionSummary;
import com.hartwig.hmftools.lilac.qc.SomaticVariantQC;
import com.hartwig.hmftools.lilac.read.BamReader;
import com.hartwig.hmftools.lilac.read.BamRecordReader;
import com.hartwig.hmftools.lilac.read.Indel;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.seq.SequenceCount;
import com.hartwig.hmftools.lilac.variant.CopyNumberAssignment;
import com.hartwig.hmftools.lilac.variant.SomaticCodingCount;
import com.hartwig.hmftools.lilac.variant.SomaticVariant;
import com.hartwig.hmftools.lilac.variant.SomaticVariantAnnotation;

public class LilacApplication
{
    private final LilacConfig mConfig;
    private final ReferenceData mRefData;

    private BamReader mRefBamReader;
    private BamReader mTumorBamReader;

    private AminoAcidFragmentPipeline mAminoAcidPipeline;
    private NucleotideFragmentFactory mNucleotideFragFactory;

    private FragmentAlleleMapper mFragAlleleMapper;
    private HlaYCoverage mHlaYCoverage;

    // key state and results
    private SequenceCount mRefAminoAcidCounts;
    private SequenceCount mRefNucleotideCounts;

    private final List<ComplexCoverage> mRankedComplexes;
    private final List<Fragment> mRefNucleotideFrags;
    private final List<FragmentAlleles> mRefFragAlleles;
    private final List<SomaticVariant> mSomaticVariants;
    private final List<SomaticCodingCount> mSomaticCodingCounts;

    private ComplexCoverage mTumorCoverage;
    private final List<Double> mTumorCopyNumber;
    private ComplexCoverage mRnaCoverage;

    private LilacQC mSummaryMetrics;
    private SolutionSummary mSolutionSummary;

    private final ResultsWriter mResultsWriter;

    public LilacApplication(final LilacConfig config, final ConfigBuilder configBuilder)
    {
        mConfig = config;
        mRefData = new ReferenceData(mConfig.ResourceDir, mConfig);

        mAminoAcidPipeline = null;
        mNucleotideFragFactory = null;
        mFragAlleleMapper = null;
        mHlaYCoverage = null;

        mRefBamReader = null;
        mTumorBamReader = null;

        // key state and results
        mRefAminoAcidCounts = null;
        mRefNucleotideCounts = null;

        mRankedComplexes = Lists.newArrayList();
        mRefNucleotideFrags = Lists.newArrayList();
        mRefFragAlleles = Lists.newArrayList();

        mSomaticVariants = Lists.newArrayList();
        mSomaticCodingCounts = Lists.newArrayList();

        mTumorCoverage = null;
        mTumorCopyNumber = Lists.newArrayList(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        mRnaCoverage = null;

        mResultsWriter = new ResultsWriter(mConfig, configBuilder);
    }

    public void setBamReaders(final BamReader refBamReader, final BamReader tumorBamReader)
    {
        mRefBamReader = refBamReader;
        mTumorBamReader = tumorBamReader;
    }

    public ReferenceData getReferenceData() { return mRefData; }

    public LilacQC getSummaryMetrics() { return mSummaryMetrics; }
    public SolutionSummary getSolutionSummary() { return mSolutionSummary; }
    public List<ComplexCoverage> getRankedComplexes() { return mRankedComplexes; }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();

        mConfig.logParams();

        if(!mRefData.load())
        {
            LL_LOGGER.error("reference data loading failed");
            System.exit(1);
        }

        boolean allValid = true;

        String referenceBam = mConfig.tumorOnly() ? mConfig.TumorBam : mConfig.ReferenceBam;

        LL_LOGGER.info("finding read support in {} bam {}", mConfig.tumorOnly() ? "tumor" : "reference", referenceBam);

        mNucleotideFragFactory = new NucleotideFragmentFactory(mRefData);

        if(mRefBamReader == null)
            mRefBamReader = new BamRecordReader(referenceBam, mConfig, GENE_CACHE.GeneTranscriptMap, mNucleotideFragFactory);

        if(mTumorBamReader == null)
        {
            if(mConfig.tumorOnly())
                mTumorBamReader = mRefBamReader;
            else if(!mConfig.TumorBam.isEmpty())
                mTumorBamReader = new BamRecordReader(mConfig.TumorBam, mConfig, GENE_CACHE.GeneTranscriptMap, mNucleotideFragFactory);
        }

        List<Fragment> refFragments = mRefBamReader.findGeneFragments();

        NUC_GENE_FRAG_ENRICHMENT.checkAddAdditionalGenes(refFragments);

        mRefNucleotideFrags.addAll(refFragments);

        byte medianBaseQuality = mNucleotideFragFactory.calculatePercentileBaseQuality(mRefNucleotideFrags,
                BASE_QUAL_PERCENTILE);

        if(medianBaseQuality < LOW_BASE_QUAL_THRESHOLD)
        {
            LL_LOGGER.info("lowering min base quality({}) to median({})", LOW_BASE_QUAL_THRESHOLD,
                    medianBaseQuality);
            LOW_BASE_QUAL_THRESHOLD = medianBaseQuality;
        }

        final Map<HlaGene, int[]> geneBaseDepth = calculateGeneCoverage(mRefNucleotideFrags);
        if(!hasSufficientGeneDepth(geneBaseDepth))
        {
            mResultsWriter.writeFailedSampleFileOutputs(geneBaseDepth, medianBaseQuality);
            return;
        }

        allValid &= validateFragments(mRefNucleotideFrags);

        mAminoAcidPipeline = new AminoAcidFragmentPipeline(mRefNucleotideFrags);

        List<Fragment> refAminoAcidFrags = mAminoAcidPipeline.highQualRefFragments();
        int totalFragmentCount = refAminoAcidFrags.size();

        LL_LOGGER.info(format("totalFrags(%d)", totalFragmentCount));

        Candidates candidateFactory = new Candidates(mRefData.NucleotideSequences, mRefData.AminoAcidSequences);

        List<GeneTask> geneTasks = Lists.newArrayList();
        geneTasks.add(
                new GeneTask(mConfig, mRefData, mAminoAcidPipeline, candidateFactory, HLA_CONTEXT_FACTORY.hlaA()));
        geneTasks.add(
                new GeneTask(mConfig, mRefData, mAminoAcidPipeline, candidateFactory, HLA_CONTEXT_FACTORY.hlaB()));
        geneTasks.add(
                new GeneTask(mConfig, mRefData, mAminoAcidPipeline, candidateFactory, HLA_CONTEXT_FACTORY.hlaC()));

        List<Callable<Void>> callableList = Lists.newArrayList(geneTasks);

        if(!TaskExecutor.executeTasks(callableList, mConfig.Threads))
            System.exit(1);

        List<HlaAllele> candidateAlleles = Lists.newArrayList();

        geneTasks.forEach(x -> x.addPhasedCandidates(candidateAlleles));

        List<PhasedEvidence> combinedPhasedEvidence = Lists.newArrayList();
        geneTasks.forEach(x -> combinedPhasedEvidence.addAll(x.phasedEvidence()));

        List<HlaAllele> recoveredAlleles = Lists.newArrayList();

        // make special note of the known stop-loss INDEL on HLA-C
        Map<HlaAllele, List<Fragment>> knownStopLossFragments = Maps.newHashMap();

        for(Map.Entry<Indel, List<Fragment>> entry : mRefBamReader.getKnownStopLossFragments().entrySet())
        {
            HlaAllele allele = mRefData.KnownStopLossIndelAlleles.get(entry.getKey());

            if(allele != null)
            {
                LL_LOGGER.info("recovering stop loss allele({}) with {} fragments", allele,
                        entry.getValue().size());

                knownStopLossFragments.put(allele, entry.getValue());

                if(!candidateAlleles.contains(allele))
                    recoveredAlleles.add(allele);
            }
        }

        if(mConfig.RestrictedAlleles.isEmpty())
        {
            // add common alleles - either if supported or forced for inclusion
            List<HlaAllele> missedCommonAlleles = mRefData.CommonAlleles.stream()
                    .filter(x -> !candidateAlleles.contains(x))
                    .filter(x -> !recoveredAlleles.contains(x))
                    .collect(Collectors.toList());

            if(!missedCommonAlleles.isEmpty())
            {
                Collections.sort(missedCommonAlleles);
                LL_LOGGER.info("recovering {} common alleles: {}", missedCommonAlleles.size(), HlaAllele.toString(missedCommonAlleles));
                recoveredAlleles.addAll(missedCommonAlleles);
            }
        }

        candidateAlleles.addAll(recoveredAlleles);

        List<HlaAllele> missingExpected = mConfig.ActualAlleles.stream()
                .filter(x -> !candidateAlleles.contains(x))
                .collect(Collectors.toList());

        if(!missingExpected.isEmpty())
        {
            LL_LOGGER.info("  including {} known actual alleles as candidates: {}",
                    missingExpected.size(), HlaAllele.toString(missingExpected));

            candidateAlleles.addAll(missingExpected);
        }

        List<HlaSequenceLoci> candidateSequences = mRefData.AminoAcidSequences.stream()
                .filter(x -> candidateAlleles.contains(x.Allele)).collect(Collectors.toList());

        // calculate allele coverage
        mRefAminoAcidCounts = SequenceCount.buildFromAminoAcids(MIN_EVIDENCE_FACTOR, refAminoAcidFrags);
        mRefNucleotideCounts = SequenceCount.buildFromNucleotides(MIN_EVIDENCE_FACTOR, refAminoAcidFrags);

        Map<HlaGene, List<Integer>> refNucleotideHetLociMap = calcNucleotideHeterogygousLoci(
                Lists.newArrayList(mRefNucleotideCounts.heterozygousLoci()));

        List<HlaSequenceLoci> candidateNucSequences = mRefData.NucleotideSequences.stream()
                .filter(x -> candidateAlleles.contains(x.Allele.asFourDigit())).collect(Collectors.toList());

        List<HlaSequenceLoci> recoveredSequences = mRefData.AminoAcidSequences.stream()
                .filter(x -> recoveredAlleles.contains(x.Allele)).collect(Collectors.toList());

        Map<HlaGene, Map<Integer, Set<String>>> geneAminoAcidHetLociMap =
                extractHeterozygousLociSequences(mAminoAcidPipeline.getReferenceAminoAcidCounts(), recoveredSequences);

        mFragAlleleMapper = new FragmentAlleleMapper(
                geneAminoAcidHetLociMap, refNucleotideHetLociMap, mAminoAcidPipeline.getReferenceNucleotides());

        mFragAlleleMapper.setKnownStopLossAlleleFragments(knownStopLossFragments);

        mRefFragAlleles.addAll(mFragAlleleMapper.createFragmentAlleles(refAminoAcidFrags, candidateSequences,
                candidateNucSequences));

        if(mRefFragAlleles.isEmpty())
        {
            LL_LOGGER.fatal("failed to assign fragments to alleles");
            System.exit(1);
        }

        mHlaYCoverage = new HlaYCoverage(mRefData.HlaYNucleotideSequences, geneAminoAcidHetLociMap, mConfig);
        mHlaYCoverage.checkThreshold(mRefFragAlleles, refAminoAcidFrags);

        // build and score complexes
        ComplexBuilder complexBuilder = new ComplexBuilder(mConfig, mRefData);

        complexBuilder.filterCandidates(mRefFragAlleles, candidateAlleles, recoveredAlleles);
        allValid &= validateAlleles(complexBuilder.getUniqueProteinAlleles());

        // reassess fragment-allele assignment with the reduced set of filtered alleles
        List<HlaAllele> confirmedRecoveredAlleles = complexBuilder.getConfirmedRecoveredAlleles();
        recoveredSequences = recoveredSequences.stream()
                .filter(x -> confirmedRecoveredAlleles.contains(x.Allele)).collect(Collectors.toList());

        geneAminoAcidHetLociMap = extractHeterozygousLociSequences(
                mAminoAcidPipeline.getReferenceAminoAcidCounts(), recoveredSequences);

        mFragAlleleMapper.setHetAminoAcidLoci(geneAminoAcidHetLociMap);
        mHlaYCoverage.updateAminoAcidLoci(geneAminoAcidHetLociMap);

        List<HlaAllele> confirmedAlleles = complexBuilder.getUniqueProteinAlleles();

        candidateSequences = candidateSequences.stream()
                .filter(x -> confirmedAlleles.contains(x.Allele)).collect(Collectors.toList());

        candidateNucSequences = candidateNucSequences.stream()
                .filter(x -> confirmedAlleles.contains(x.Allele.asFourDigit())).collect(Collectors.toList());

        mRefFragAlleles.clear();

        LL_LOGGER.debug(
                "creating fragment alleles from amnioAcidfrags({}) candidateSequences({}) nucleotideSequences({})",
                refAminoAcidFrags.size(), candidateSequences.size(), candidateNucSequences.size());

        mRefFragAlleles.addAll(mFragAlleleMapper.createFragmentAlleles(refAminoAcidFrags, candidateSequences,
                candidateNucSequences));

        // down-sample if ref depth is higher than configured cap
        List<FragmentAlleles> calcRefFragAlleles = checkDownsampleRefFragmentAlleles();

        List<HlaComplex> complexes = complexBuilder.buildComplexes(calcRefFragAlleles, confirmedRecoveredAlleles);
        // allValid &= validateComplexes(complexes); // too expensive in current form even for validation, address in unit tests instead

        ComplexCoverageCalculator complexCalculator = new ComplexCoverageCalculator(mConfig);

        LL_LOGGER.info("calculating coverage for complexes({}) and ref alleles({})",
                complexes.size(),
                mRefFragAlleles.size() > calcRefFragAlleles.size()
                        ? format("%d capped=%d", mRefFragAlleles.size(), calcRefFragAlleles.size())
                        : mRefFragAlleles.size());

        List<ComplexCoverage> calculatedComplexes = complexCalculator.calculateComplexCoverages(calcRefFragAlleles, complexes);

        ComplexCoverageRanking complexRanker = new ComplexCoverageRanking(mConfig.TopScoreThreshold, mRefData);
        mRankedComplexes.addAll(complexRanker.rankCandidates(calculatedComplexes, recoveredAlleles, candidateSequences));

        if(mRankedComplexes.isEmpty())
        {
            LL_LOGGER.fatal("failed to calculate complex coverage");
            System.exit(1);
        }

        if(calcRefFragAlleles.size() < mRefFragAlleles.size())
        {
            List<HlaComplex> filteredComplexes = mRankedComplexes.stream().map(ComplexCoverage::toComplex).toList();

            LL_LOGGER.debug("recalculating coverage for complexes({}) and ref alleles({})",
                    filteredComplexes.size(), mRefFragAlleles.size());

            calculatedComplexes = complexCalculator.calculateComplexCoverages(mRefFragAlleles, filteredComplexes);
            complexRanker = new ComplexCoverageRanking(0, mRefData);
            mRankedComplexes.clear();
            mRankedComplexes.addAll(complexRanker.rankCandidates(calculatedComplexes, recoveredAlleles, candidateSequences));
        }

        ComplexCoverage winningRefCoverage = mRankedComplexes.get(0);

        if(!mConfig.ActualAlleles.isEmpty())
        {
            // search amongst the ranked candidates otherwise manually compute coverage
            ComplexCoverage actualAllelesCoverage = null;

            for(ComplexCoverage rankedCoverage : mRankedComplexes)
            {
                if(rankedCoverage.getAlleles().size() != mConfig.ActualAlleles.size())
                    continue;

                if(mConfig.ActualAlleles.containsAll(rankedCoverage.getAlleles()))
                {
                    actualAllelesCoverage = rankedCoverage;
                    break;
                }
            }

            if(actualAllelesCoverage == null)
            {
                actualAllelesCoverage = ComplexBuilder.calcProteinCoverage(mRefFragAlleles, mConfig.ActualAlleles);
            }

            LL_LOGGER.info("forcing winning coverage to actual alleles");
            winningRefCoverage = actualAllelesCoverage;
        }

        winningRefCoverage.expandToSixAlleles();

        List<HlaAllele> winningAlleles = winningRefCoverage.getAlleles();

        mHlaYCoverage.assignReferenceFragments(winningAlleles, mRefFragAlleles, refAminoAcidFrags);

        List<HlaSequenceLoci> winningSequences = candidateSequences.stream()
                .filter(x -> winningAlleles.contains(x.Allele)).collect(Collectors.toList());

        List<HlaSequenceLoci> winningNucSequences = candidateNucSequences.stream()
                .filter(x -> winningAlleles.contains(x.Allele.asFourDigit()))
                .collect(Collectors.toList());

        LL_LOGGER.info("{}", HlaComplexFile.header());

        for(ComplexCoverage rankedComplex : mRankedComplexes)
        {
            LL_LOGGER.info(HlaComplexFile.asString(rankedComplex));
        }

        // log key results for fast post-run analysis
        StringJoiner totalCoverages = new StringJoiner(",");
        winningRefCoverage.getAlleleCoverage().forEach(x -> totalCoverages.add(format("%.0f", x.TotalCoverage)));

        double scoreMargin = 0;
        StringJoiner nextSolutionInfo = new StringJoiner(ITEM_DELIM);

        if(mRankedComplexes.size() > 1)
        {
            ComplexCoverage nextSolution = mRankedComplexes.get(1);
            scoreMargin = mRankedComplexes.get(0).getScore() - nextSolution.getScore();
            nextSolution.getAlleles().stream().filter(x -> !winningAlleles.contains(x))
                    .forEach(x -> nextSolutionInfo.add(x.toString()));
        }

        LL_LOGGER.info("WINNERS_REF: {}, {}, {}, {}, {}, {}",
                mConfig.Sample, mRankedComplexes.size(), HlaAllele.toString(winningRefCoverage.getAlleles()), format("%.3f", scoreMargin),
                nextSolutionInfo, totalCoverages);

        // write fragment assignment data
        for(FragmentAlleles fragAllele : mRefFragAlleles)
        {
            if(fragAllele.getFragment().isScopeSet())
                continue;

            if(winningAlleles.stream().anyMatch(fragAllele::contains))
                fragAllele.getFragment().setScope(SOLUTION);
            else
                fragAllele.getFragment().setScope(CANDIDATE);
        }

        mSomaticCodingCounts.addAll(SomaticCodingCount.create(winningAlleles));

        extractTumorResults(winningAlleles, winningRefCoverage, winningSequences, winningNucSequences);

        extractRnaCoverage(winningAlleles, winningSequences, winningNucSequences);

        if(!allValid)
        {
            LL_LOGGER.error("failed validation");
            mResultsWriter.writeFailedSampleFileOutputs(geneBaseDepth, medianBaseQuality);
            return;
        }

        // create various QC and other metrics
        SomaticVariantQC somaticVariantQC = SomaticVariantQC.create(mSomaticVariants.size(), mSomaticCodingCounts);

        List<Fragment> unmatchedFrags = refAminoAcidFrags.stream().filter(x -> x.scope().isUnmatched())
                .collect(Collectors.toList());

        HaplotypeQC haplotypeQC = HaplotypeQC.create(
                winningSequences, mRefData.HlaYAminoAcidSequences, combinedPhasedEvidence, mRefAminoAcidCounts,
                unmatchedFrags);

        AminoAcidQC aminoAcidQC = AminoAcidQC.create(
                winningSequences, mRefData.HlaYAminoAcidSequences, mRefAminoAcidCounts,
                haplotypeQC.UnmatchedHaplotypes, totalFragmentCount);

        BamQC bamQC = BamQC.create(mRefBamReader, geneBaseDepth);
        CoverageQC coverageQC = CoverageQC.create(refAminoAcidFrags, winningRefCoverage);

        mSummaryMetrics = new LilacQC(
                scoreMargin, nextSolutionInfo.toString(), medianBaseQuality, mHlaYCoverage.getSelectedAllele(),
                aminoAcidQC, bamQC, coverageQC, haplotypeQC, somaticVariantQC);

        ComplexCoverage refCoverage = !mConfig.tumorOnly() ? winningRefCoverage
                : ComplexCoverage.create(Lists.newArrayList());

        mSolutionSummary = SolutionSummary.create(refCoverage, mTumorCoverage, mTumorCopyNumber,
                mSomaticCodingCounts, mRnaCoverage);

        writeFileOutputs();

        mResultsWriter.close();

        LL_LOGGER.info("Lilac complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private List<FragmentAlleles> checkDownsampleRefFragmentAlleles()
    {
        if(mConfig.MaxRefFragments == 0 || mRefFragAlleles.size() <= mConfig.MaxRefFragments)
            return mRefFragAlleles;

        if(mRefFragAlleles.size() <= mConfig.MaxRefFragments * 2)
        {
            return mRefFragAlleles.stream().sorted(Comparator.comparing(x -> x.getFragment().id()))
                    .limit(mConfig.MaxRefFragments)
                    .collect(Collectors.toList());
        }

        List<FragmentAlleles> calcRefFragAlleles = Lists.newArrayList();
        int nthElement = (int) floor(mRefFragAlleles.size() / (double) mConfig.MaxRefFragments);

        int counter = 0;
        for(FragmentAlleles fragmentAllele : mRefFragAlleles)
        {
            ++counter;

            if(counter == nthElement)
            {
                calcRefFragAlleles.add(fragmentAllele);

                if(calcRefFragAlleles.size() >= mConfig.MaxRefFragments)
                    break;

                counter = 0;
            }
        }

        return calcRefFragAlleles;
    }

    private void extractTumorResults(
            final List<HlaAllele> winningAlleles, final ComplexCoverage winningRefCoverage,
            final List<HlaSequenceLoci> winningSequences, final List<HlaSequenceLoci> winningNucSequences)
    {
        if(mConfig.TumorBam.isEmpty())
        {
            mTumorCoverage = ComplexCoverage.create(Lists.newArrayList());
            return;
        }

        if(mConfig.tumorOnly())
        {
            mTumorCoverage = winningRefCoverage;
        }
        else
        {
            List<Fragment> tumorNucleotideFrags = mTumorBamReader.findGeneFragments();
            NUC_GENE_FRAG_ENRICHMENT.checkAddAdditionalGenes(tumorNucleotideFrags);

            List<Fragment> tumorFragments = mAminoAcidPipeline.calcComparisonCoverageFragments(tumorNucleotideFrags);

            mResultsWriter.writeFragments(TUMOR, tumorFragments);

            LL_LOGGER.info("calculating tumor coverage from frags({} highQual={})", tumorNucleotideFrags.size(), tumorFragments.size());

            List<FragmentAlleles> tumorFragAlleles =
                    mFragAlleleMapper.createFragmentAlleles(tumorFragments, winningSequences, winningNucSequences);

            if(mHlaYCoverage.exceedsThreshold())
                mHlaYCoverage.assignFragments(winningAlleles, tumorFragAlleles, tumorFragments, TUMOR);

            mTumorCoverage = ComplexBuilder.calcProteinCoverage(tumorFragAlleles, winningAlleles);

            mTumorCoverage.populateMissingCoverage(winningAlleles);
        }

        if(!mConfig.CopyNumberFile.isEmpty())
        {
            CopyNumberAssignment copyNumberAssignment = new CopyNumberAssignment();
            copyNumberAssignment.loadCopyNumberData(mConfig);

            copyNumberAssignment.assign(mConfig.Sample, winningRefCoverage, mTumorCoverage, mTumorCopyNumber);
        }

        // SOMATIC VARIANTS
        SomaticVariantAnnotation variantAnnotation = new SomaticVariantAnnotation(mConfig, GENE_CACHE.GeneTranscriptMap);

        if(variantAnnotation.getSomaticVariants().isEmpty())
            return;

        mSomaticVariants.addAll(variantAnnotation.getSomaticVariants());

        LL_LOGGER.info("calculating somatic variant allele coverage");

        for(SomaticVariant variant : mSomaticVariants)
        {
            List<Fragment> variantFragments = mTumorBamReader.findVariantFragments(variant);

            List<AlleleCoverage> variantCoverage = variantAnnotation.assignAlleleCoverage(variant, variantFragments, winningSequences);

            List<HlaAllele> variantAlleles = variantCoverage.stream().map(x -> x.Allele).collect(Collectors.toList());
            LL_LOGGER.info("  {} -> {}}", variant, variantCoverage);

            mResultsWriter.writeVariant(variant.Context, variantAlleles);

            addVariant(mSomaticCodingCounts, variant, variantAlleles);
        }
    }

    private void extractRnaCoverage(final List<HlaAllele> winningAlleles, final List<HlaSequenceLoci> winningSequences,
            final List<HlaSequenceLoci> winningNucSequences)
    {
        mRnaCoverage = LilacAppendRna.extractRnaCoverage(
                mConfig.RnaBam, mConfig, mRefData, mNucleotideFragFactory, NUC_GENE_FRAG_ENRICHMENT, mAminoAcidPipeline, mFragAlleleMapper,
                winningAlleles, winningSequences, winningNucSequences);
    }

    private void writeFileOutputs()
    {
        mSummaryMetrics.log(mConfig.Sample);

        LL_LOGGER.info("writing output to {}", mConfig.OutputDir);

        mResultsWriter.writeMainOutputs(mSummaryMetrics, mSolutionSummary, mRankedComplexes);

        mResultsWriter.writeDetailedOutputs(mRefAminoAcidCounts, mRefNucleotideCounts, mAminoAcidPipeline, mHlaYCoverage);

        mResultsWriter.writeReferenceFragments(mRankedComplexes, mRefNucleotideFrags, mRefFragAlleles);
    }

    private boolean validateFragments(final Collection<Fragment> fragments)
    {
        if(!mConfig.RunValidation)
            return true;

        List<Fragment> invalidFragments = fragments.stream().filter(x -> !x.validate()).toList();
        if(invalidFragments.isEmpty())
            return true;

        LL_LOGGER.warn("has {} invalid fragments", invalidFragments.size());
        return false;
    }

    private boolean validateAlleles(final List<HlaAllele> alleles)
    {
        // check uniqueness amongst valid candidate alleles
        if(!mConfig.RunValidation)
            return true;

        Set<HlaAllele> duplicateAlleles = HlaAllele.findDuplicates(alleles);
        if(duplicateAlleles.isEmpty())
            return true;

        LL_LOGGER.warn("has {} duplicate alleles from complex building", duplicateAlleles.size());
        return false;
    }

    private boolean validateComplexes(final List<HlaComplex> complexes)
    {
        if(!mConfig.RunValidation)
            return true;

        Set<HlaComplex> duplicates = findDuplicates(complexes);
        if(duplicates.isEmpty())
            return true;

        LL_LOGGER.warn("has {} duplicate complexes", duplicates.size());
        return false;
    }

    private boolean hasSufficientGeneDepth(final Map<HlaGene, int[]> geneBaseDepth)
    {
        int aLowCoveragePositions = (int) Arrays.stream(geneBaseDepth.get(HLA_A)).filter(x -> x < WARN_LOW_COVERAGE_DEPTH).count();
        int bLowCoveragePositions = (int) Arrays.stream(geneBaseDepth.get(HLA_B)).filter(x -> x < WARN_LOW_COVERAGE_DEPTH).count();
        int cLowCoveragePositions = (int) Arrays.stream(geneBaseDepth.get(HLA_C)).filter(x -> x < WARN_LOW_COVERAGE_DEPTH).count();

        int totalLowCoveragePositions = aLowCoveragePositions + bLowCoveragePositions + cLowCoveragePositions;

        if(!mConfig.ReferenceBam.isEmpty() && totalLowCoveragePositions >= mConfig.FatalTotalLowCoveragePositions)
        {
            LL_LOGGER.warn("exiting due to too low coverage: bases with <{} coverage per gene(A={}, B={}, C={}) total({}) exceeds threshold({})",
                    WARN_LOW_COVERAGE_DEPTH, aLowCoveragePositions, bLowCoveragePositions, cLowCoveragePositions,
                    totalLowCoveragePositions, mConfig.FatalTotalLowCoveragePositions);
            return false;
        }

        return true;
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        LilacConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        LilacConfig lilacConfig = new LilacConfig(configBuilder);
        LilacApplication lilac = new LilacApplication(lilacConfig, configBuilder);

        File stackSampleFile = lilacConfig.OutputDir.isEmpty()
                ? null
                : new File(new File(lilacConfig.OutputDir), format("%s.lilac.stacks", lilacConfig.Sample));
        try(StackSampler stackSampler = stackSampleFile == null || lilacConfig.StackSampleRate <= 0
                ? null
                : new StackSampler(lilacConfig.StackSampleRate, stackSampleFile))
        {
            lilac.run();
        }
        catch(Exception e)
        {
            throw new RuntimeException(e);
        }
    }
}
