package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConfig.LOG_DEBUG;
import static com.hartwig.hmftools.lilac.LilacConfig.LOG_LEVEL;
import static com.hartwig.hmftools.lilac.LilacConstants.A_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.LilacConstants.B_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.LilacConstants.C_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_A;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_B;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_C;
import static com.hartwig.hmftools.lilac.LilacConstants.ITEM_DELIM;
import static com.hartwig.hmftools.lilac.evidence.Candidates.addPhasedCandidates;
import static com.hartwig.hmftools.lilac.seq.SequenceCount.extractHeterozygousLociSequences;
import static com.hartwig.hmftools.lilac.evidence.NucleotideFiltering.calcNucleotideHeterogygousLoci;
import static com.hartwig.hmftools.lilac.coverage.HlaComplex.findDuplicates;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.CANDIDATE;
import static com.hartwig.hmftools.lilac.fragment.FragmentScope.SOLUTION;
import static com.hartwig.hmftools.lilac.variant.SomaticCodingCount.addVariant;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.lilac.evidence.Candidates;
import com.hartwig.hmftools.lilac.coverage.FragmentAlleleMapper;
import com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage;
import com.hartwig.hmftools.lilac.coverage.HlaComplex;
import com.hartwig.hmftools.lilac.coverage.HlaComplexBuilder;
import com.hartwig.hmftools.lilac.coverage.HlaComplexCoverage;
import com.hartwig.hmftools.lilac.coverage.ComplexCoverageCalculator;
import com.hartwig.hmftools.lilac.coverage.HlaComplexCoverageRanking;
import com.hartwig.hmftools.lilac.coverage.HlaComplexFile;
import com.hartwig.hmftools.lilac.evidence.PhasedEvidence;
import com.hartwig.hmftools.lilac.evidence.PhasedEvidenceFactory;
import com.hartwig.hmftools.lilac.fragment.AminoAcidFragmentPipeline;
import com.hartwig.hmftools.lilac.fragment.FragmentUtils;
import com.hartwig.hmftools.lilac.fragment.NucleotideGeneEnrichment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaContext;
import com.hartwig.hmftools.lilac.hla.HlaContextFactory;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.fragment.NucleotideFragmentFactory;
import com.hartwig.hmftools.lilac.read.BamReader;
import com.hartwig.hmftools.lilac.read.Indel;
import com.hartwig.hmftools.lilac.seq.SequenceCount;
import com.hartwig.hmftools.lilac.variant.CopyNumberAssignment;
import com.hartwig.hmftools.lilac.qc.AminoAcidQC;
import com.hartwig.hmftools.lilac.qc.BamQC;
import com.hartwig.hmftools.lilac.qc.CoverageQC;
import com.hartwig.hmftools.lilac.qc.HaplotypeQC;
import com.hartwig.hmftools.lilac.qc.SolutionSummary;
import com.hartwig.hmftools.lilac.qc.LilacQC;
import com.hartwig.hmftools.lilac.qc.SomaticVariantQC;
import com.hartwig.hmftools.lilac.coverage.FragmentAlleles;
import com.hartwig.hmftools.lilac.read.BamRecordReader;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.variant.LilacVCF;
import com.hartwig.hmftools.lilac.variant.SomaticCodingCount;
import com.hartwig.hmftools.lilac.variant.SomaticVariant;
import com.hartwig.hmftools.lilac.variant.SomaticVariantAnnotation;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class LilacApplication
{
    private final LilacConfig mConfig;
    private final ReferenceData mRefData;

    private BamReader mRefBamReader;
    private BamReader mTumorBamReader;

    private AminoAcidFragmentPipeline mAminoAcidPipeline;
    private final NucleotideGeneEnrichment mNucleotideGeneEnrichment;
    private NucleotideFragmentFactory mNucleotideFragFactory;

    // key state and results
    private SequenceCount mRefAminoAcidCounts;
    private SequenceCount mRefNucleotideCounts;

    private List<HlaComplexCoverage> mRankedComplexes;
    private List<Fragment> mRefNucleotideFrags;
    private List<FragmentAlleles> mRefFragAlleles;
    private final List<SomaticVariant> mSomaticVariants;
    private final List<SomaticCodingCount> mSomaticCodingCounts;

    private HlaComplexCoverage mTumorCoverage;
    private final List<Double> mTumorCopyNumber;
    private HlaComplexCoverage mRnaCoverage;

    private LilacQC mSummaryMetrics;
    private SolutionSummary mSolutionSummary;

    public LilacApplication(final LilacConfig config)
    {
        mConfig = config;
        mRefData = new ReferenceData(mConfig.ResourceDir, mConfig);

        mAminoAcidPipeline = null;
        mNucleotideGeneEnrichment = new NucleotideGeneEnrichment(A_EXON_BOUNDARIES, B_EXON_BOUNDARIES, C_EXON_BOUNDARIES);
        mNucleotideFragFactory = null;

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
    }

    public void setBamReaders(final BamReader refBamReader, final BamReader tumorBamReader)
    {
        mRefBamReader = refBamReader;
        mTumorBamReader = tumorBamReader;
    }

    public ReferenceData getReferenceData() { return mRefData; }

    public LilacQC getSummaryMetrics() { return mSummaryMetrics; }
    public SolutionSummary getSolutionSummary() { return mSolutionSummary; }
    public List<HlaComplexCoverage> getRankedComplexes() { return mRankedComplexes; }

    public void run()
    {
        if(!mConfig.isValid())
        {
            LL_LOGGER.warn("invalid config, exiting");
            System.exit(1);
        }

        final VersionInfo version = new VersionInfo("lilac.version");
        LL_LOGGER.info("Lilac version({}), key parameters:", version.version());
        mConfig.logParams();

        HlaContextFactory hlaContextFactory = new HlaContextFactory(A_EXON_BOUNDARIES, B_EXON_BOUNDARIES, C_EXON_BOUNDARIES);
        HlaContext hlaAContext = hlaContextFactory.hlaA();
        HlaContext hlaBContext = hlaContextFactory.hlaB();
        HlaContext hlaCContext = hlaContextFactory.hlaC();

        if(!mRefData.load())
        {
            LL_LOGGER.error("reference data loading failed");
            System.exit(1);
        }

        boolean allValid = true;

        LL_LOGGER.info("querying records from reference bam {}", mConfig.ReferenceBam);

        mNucleotideFragFactory = new NucleotideFragmentFactory(
                mConfig.MinBaseQual, mRefData.AminoAcidSequencesWithInserts, mRefData.AminoAcidSequencesWithDeletes,
                mRefData.LociPositionFinder);

        if(mRefBamReader == null)
            mRefBamReader = new BamRecordReader(mConfig.ReferenceBam, mConfig.RefGenome, mRefData.HlaTranscriptData, mNucleotideFragFactory);

        if(mTumorBamReader == null)
            mTumorBamReader = new BamRecordReader(mConfig.TumorBam, mConfig.RefGenome, mRefData.HlaTranscriptData, mNucleotideFragFactory);

        mRefNucleotideFrags.addAll(mNucleotideGeneEnrichment.enrich(mRefBamReader.readFromBam()));

        allValid &= validateFragments(mRefNucleotideFrags);

        mAminoAcidPipeline = new AminoAcidFragmentPipeline(mConfig, mRefNucleotideFrags);

        List<Fragment> refAminoAcidFrags = mAminoAcidPipeline.getReferenceFragments();
        int totalFragmentCount = refAminoAcidFrags.size();

        double minEvidence = mConfig.calcMinEvidence(totalFragmentCount);

        LL_LOGGER.info(String.format("totalFrags(%d) minEvidence(%.1f) minHighQualEvidence(%.1f)",
                totalFragmentCount, minEvidence, mConfig.calcMinHighQualEvidence(totalFragmentCount)));

        // apply special filtering and splice checks on fragments, just for use in phasing
        List<Fragment> aCandidateFrags = mAminoAcidPipeline.referencePhasingFragments(hlaAContext);
        List<Fragment> bCandidateFrags = mAminoAcidPipeline.referencePhasingFragments(hlaBContext);
        List<Fragment> cCandidateFrags = mAminoAcidPipeline.referencePhasingFragments(hlaCContext);

        // determine un-phased Candidates
        Candidates candidateFactory = new Candidates(minEvidence, mRefData.NucleotideSequences, mRefData.AminoAcidSequences);
        List<HlaAllele> aUnphasedCandidates = candidateFactory.unphasedCandidates(hlaAContext, aCandidateFrags);
        List<HlaAllele> bUnphasedCandidates = candidateFactory.unphasedCandidates(hlaBContext, bCandidateFrags);
        List<HlaAllele> cUnphasedCandidates = candidateFactory.unphasedCandidates(hlaCContext, cCandidateFrags);

        // determine phasing of amino acids
        PhasedEvidenceFactory phasedEvidenceFactory = new PhasedEvidenceFactory(mConfig, minEvidence);
        List<PhasedEvidence> aPhasedEvidence = phasedEvidenceFactory.evidence(hlaAContext, aCandidateFrags);
        List<PhasedEvidence> bPhasedEvidence = phasedEvidenceFactory.evidence(hlaBContext, bCandidateFrags);
        List<PhasedEvidence> cPhasedEvidence = phasedEvidenceFactory.evidence(hlaCContext, cCandidateFrags);

        // validate phasing against expected sequences
        List<HlaSequenceLoci> expectedSequences = Lists.newArrayList();
        if(!mConfig.ExpectedAlleles.isEmpty())
        {
            expectedSequences.addAll(mRefData.AminoAcidSequences.stream()
                    .filter(x -> mConfig.ExpectedAlleles.contains(x.Allele.asFourDigit())).collect(Collectors.toList()));

            PhasedEvidence.logInconsistentEvidence(GENE_A, aPhasedEvidence, expectedSequences);
            PhasedEvidence.logInconsistentEvidence(GENE_B, bPhasedEvidence, expectedSequences);
            PhasedEvidence.logInconsistentEvidence(GENE_C, cPhasedEvidence, expectedSequences);
        }

        // gather up all phased candidates
        List<HlaAllele> candidateAlleles = Lists.newArrayList();
        List<HlaAllele> aCandidates = candidateFactory.phasedCandidates(hlaAContext, aUnphasedCandidates, aPhasedEvidence);
        List<HlaAllele> bCandidates = candidateFactory.phasedCandidates(hlaBContext, bUnphasedCandidates, bPhasedEvidence);
        List<HlaAllele> cCandidates = candidateFactory.phasedCandidates(hlaCContext, cUnphasedCandidates, cPhasedEvidence);

        addPhasedCandidates(candidateAlleles, aCandidates, mConfig, mRefData);
        addPhasedCandidates(candidateAlleles, bCandidates, mConfig, mRefData);
        addPhasedCandidates(candidateAlleles, cCandidates, mConfig, mRefData);

        List<HlaAllele> recoveredAlleles = Lists.newArrayList();

        // make special note of the known stop-loss INDEL on HLA-C
        Map<HlaAllele,List<Fragment>> knownStopLossFragments = Maps.newHashMap();

        for(Map.Entry<Indel,List<Fragment>> entry : mRefBamReader.getKnownStopLossFragments().entrySet())
        {
            HlaAllele allele = mRefData.KnownStopLossIndelAlleles.get(entry.getKey());

            if(allele != null)
            {
                LL_LOGGER.info("recovering stop loss allele({}) with {} fragments", allele, entry.getValue().size());

                knownStopLossFragments.put(allele, entry.getValue());

                if(!candidateAlleles.contains(allele))
                    recoveredAlleles.add(allele);
            }
        }

        if(!mConfig.RestrictedAlleles.isEmpty())
        {
            // if using a set of restricted allels, make no concession for any missed or common alleles
            mConfig.ExpectedAlleles.stream().filter(x -> !candidateAlleles.contains(x)).forEach(x -> candidateAlleles.add(x));
        }
        else
        {
            // add common alleles - either if supported or forced for inclusion
            List<HlaAllele> missedCommonAlleles = mRefData.CommonAlleles.stream()
                    .filter(x -> !candidateAlleles.contains(x))
                    .filter(x -> !recoveredAlleles.contains(x))
                    .collect(Collectors.toList());

            if(!missedCommonAlleles.isEmpty())
            {
                Collections.sort(missedCommonAlleles);
                LL_LOGGER.info("recovering common alleles: {}", HlaAllele.toString(missedCommonAlleles));
                recoveredAlleles.addAll(missedCommonAlleles);
            }
        }

        candidateAlleles.addAll(recoveredAlleles);

        List<HlaAllele> missingExpected = mConfig.ExpectedAlleles.stream().filter(x -> !candidateAlleles.contains(x)).collect(Collectors.toList());

        if (!missingExpected.isEmpty())
        {
            LL_LOGGER.info("  including {} expected alleles as candidates: {}",
                    missingExpected.size(), HlaAllele.toString(missingExpected));

            candidateAlleles.addAll(missingExpected);
        }

        List<HlaSequenceLoci> candidateSequences = mRefData.AminoAcidSequences.stream()
                .filter(x -> candidateAlleles.contains(x.Allele)).collect(Collectors.toList());

        // calculate allele coverage
        mRefAminoAcidCounts = SequenceCount.aminoAcids(minEvidence, refAminoAcidFrags);
        mRefNucleotideCounts = SequenceCount.nucleotides(minEvidence, refAminoAcidFrags);

        Map<String,List<Integer>> refNucleotideHetLociMap = calcNucleotideHeterogygousLoci(mRefNucleotideCounts.heterozygousLoci());

        List<HlaSequenceLoci> candidateNucSequences = mRefData.NucleotideSequences.stream()
                .filter(x -> candidateAlleles.contains(x.Allele.asFourDigit())).collect(Collectors.toList());

        List<HlaSequenceLoci> recoveredSequences = mRefData.AminoAcidSequences.stream()
                .filter(x -> recoveredAlleles.contains(x.Allele)).collect(Collectors.toList());

        Map<String,Map<Integer,Set<String>>> geneAminoAcidHetLociMap =
                extractHeterozygousLociSequences(mAminoAcidPipeline.getReferenceAminoAcidCounts(), minEvidence, recoveredSequences);

        FragmentAlleleMapper fragAlleleMapper = new FragmentAlleleMapper(
                geneAminoAcidHetLociMap, refNucleotideHetLociMap, mAminoAcidPipeline.getReferenceNucleotides());

        fragAlleleMapper.setKnownStopLossAlleleFragments(knownStopLossFragments);

        mRefFragAlleles.addAll(fragAlleleMapper.createFragmentAlleles(refAminoAcidFrags, candidateSequences, candidateNucSequences));

        if(mRefFragAlleles.isEmpty())
        {
            LL_LOGGER.fatal("failed to assign fragments to alleles");
            System.exit(1);
        }

        boolean hasHlaY = fragAlleleMapper.checkHlaYSupport(
                mRefData.HlaYNucleotideSequences, mRefFragAlleles, refAminoAcidFrags, true, "ref");

        // build and score complexes
        HlaComplexBuilder complexBuilder = new HlaComplexBuilder(mConfig, mRefData);

        complexBuilder.filterCandidates(mRefFragAlleles, candidateAlleles, recoveredAlleles);
        allValid &= validateAlleles(complexBuilder.getUniqueProteinAlleles());

        // reassess fragment-allele assignment with the reduced set of filtered alleles
        List<HlaAllele> confirmedRecoveredAlleles = complexBuilder.getConfirmedRecoveredAlleles();
        recoveredSequences = recoveredSequences.stream()
                .filter(x -> confirmedRecoveredAlleles.contains(x.Allele)).collect(Collectors.toList());

        geneAminoAcidHetLociMap =
                extractHeterozygousLociSequences(mAminoAcidPipeline.getReferenceAminoAcidCounts(), minEvidence, recoveredSequences);

        fragAlleleMapper.setHetAminoAcidLoci(geneAminoAcidHetLociMap);

        List<HlaAllele> confirmedAlleles = complexBuilder.getUniqueProteinAlleles();

        candidateSequences = candidateSequences.stream()
                .filter(x -> confirmedAlleles.contains(x.Allele)).collect(Collectors.toList());

        candidateNucSequences = candidateNucSequences.stream()
                .filter(x -> confirmedAlleles.contains(x.Allele.asFourDigit())).collect(Collectors.toList());

        mRefFragAlleles.clear();
        mRefFragAlleles.addAll(fragAlleleMapper.createFragmentAlleles(refAminoAcidFrags, candidateSequences, candidateNucSequences));

        List<HlaComplex> complexes = complexBuilder.buildComplexes(mRefFragAlleles, confirmedRecoveredAlleles);
        // allValid &= validateComplexes(complexes); // too expensive in current form even for validation, address in unit tests instead

        LL_LOGGER.info("calculating coverage of {} complexes", complexes.size());
        ComplexCoverageCalculator complexCalculator = new ComplexCoverageCalculator(mConfig.Threads);
        List<HlaComplexCoverage> calculatedComplexes = complexCalculator.calculateComplexCoverages(mRefFragAlleles, complexes);

        HlaComplexCoverageRanking complexRanker = new HlaComplexCoverageRanking(mConfig.TopScoreThreshold, mRefData);
        mRankedComplexes.addAll(complexRanker.rankCandidates(calculatedComplexes, recoveredAlleles, candidateSequences));

        if(mRankedComplexes.isEmpty())
        {
            LL_LOGGER.fatal("failed to calculate complex coverage");
            System.exit(1);
        }

        if(!expectedSequences.isEmpty())
        {
            HlaComplexCoverage expectedCoverage = HlaComplexBuilder.calcProteinCoverage(
                    mRefFragAlleles, expectedSequences.stream().map(x -> x.Allele).collect(Collectors.toList()));
            LL_LOGGER.info("expected allele coverage: {}", expectedCoverage);
        }

        HlaComplexCoverage winningRefCoverage = mRankedComplexes.get(0);
        winningRefCoverage.expandToSixAlleles();

        List<HlaAllele> winningAlleles = winningRefCoverage.getAlleles();

        List<HlaSequenceLoci> winningSequences = candidateSequences.stream()
                .filter(x -> winningAlleles.contains(x.Allele)).collect(Collectors.toList());

        List<HlaSequenceLoci> winningNucSequences = candidateNucSequences.stream()
                .filter(x -> winningAlleles.contains(x.Allele.asFourDigit())).collect(Collectors.toList());

        LL_LOGGER.info("{}", HlaComplexFile.header());

        for (HlaComplexCoverage rankedComplex : mRankedComplexes)
        {
            LL_LOGGER.info(HlaComplexFile.asString(rankedComplex));
        }

        // log key results for fast post-run analysis
        StringJoiner totalCoverages = new StringJoiner(",");
        winningRefCoverage.getAlleleCoverage().forEach(x -> totalCoverages.add(String.format("%.0f",x.TotalCoverage)));

        double scoreMargin = 0;
        StringJoiner nextSolutionInfo = new StringJoiner(ITEM_DELIM);

        if(mRankedComplexes.size() > 1)
        {
            HlaComplexCoverage nextSolution = mRankedComplexes.get(1);
            scoreMargin = mRankedComplexes.get(0).getScore() - nextSolution.getScore();
            nextSolution.getAlleles().stream().filter(x -> !winningAlleles.contains(x)).forEach(x -> nextSolutionInfo.add(x.toString()));
        }

        LL_LOGGER.info("WINNERS_{}: {}, {}, {}, {}, {}, {}",
                mConfig.RunId.isEmpty() ? "REF" : mConfig.RunId, mConfig.Sample, mRankedComplexes.size(),
                HlaAllele.toString(winningRefCoverage.getAlleles()), String.format("%.3f", scoreMargin), nextSolutionInfo, totalCoverages);

        // write fragment assignment data
        for(FragmentAlleles fragAllele : mRefFragAlleles)
        {
            if(fragAllele.getFragment().isScopeSet())
                continue;

            if(winningAlleles.stream().anyMatch(x -> fragAllele.contains(x)))
                fragAllele.getFragment().setScope(SOLUTION);
            else
                fragAllele.getFragment().setScope(CANDIDATE);
        }

        mSomaticCodingCounts.addAll(SomaticCodingCount.create(winningAlleles));

        extractTumorResults(fragAlleleMapper, winningAlleles, winningRefCoverage, winningSequences, winningNucSequences);

        extractRnaCoverage(fragAlleleMapper, winningAlleles, winningSequences, winningNucSequences);

        if(!allValid)
        {
            LL_LOGGER.error("failed validation");
        }

        // create various QC and other metrics
        SomaticVariantQC somaticVariantQC = SomaticVariantQC.create(mSomaticVariants.size(), mSomaticCodingCounts);

        List<PhasedEvidence> combinedPhasedEvidence = Lists.newArrayList();
        combinedPhasedEvidence.addAll(aPhasedEvidence);
        combinedPhasedEvidence.addAll(bPhasedEvidence);
        combinedPhasedEvidence.addAll(cPhasedEvidence);

        List<Fragment> unmatchedFrags = refAminoAcidFrags.stream().filter(x -> x.scope().isUnmatched()).collect(Collectors.toList());

        HaplotypeQC haplotypeQC = HaplotypeQC.create(
                winningSequences, mRefData.HlaYAminoAcidSequences, combinedPhasedEvidence, mRefAminoAcidCounts, unmatchedFrags);

        AminoAcidQC aminoAcidQC = AminoAcidQC.create(
                winningSequences, mRefData.HlaYAminoAcidSequences, mRefAminoAcidCounts,
                haplotypeQC.UnmatchedHaplotypes, totalFragmentCount);

        BamQC bamQC = BamQC.create(mRefBamReader);
        CoverageQC coverageQC = CoverageQC.create(refAminoAcidFrags, winningRefCoverage);

        mSummaryMetrics = new LilacQC(
                hasHlaY, scoreMargin, nextSolutionInfo.toString(), aminoAcidQC, bamQC, coverageQC, haplotypeQC, somaticVariantQC);

        mSolutionSummary = SolutionSummary.create(winningRefCoverage, mTumorCoverage, mTumorCopyNumber, mSomaticCodingCounts, mRnaCoverage);
    }

    public void extractTumorResults(
            final FragmentAlleleMapper fragAlleleMapper, final List<HlaAllele> winningAlleles, final HlaComplexCoverage winningRefCoverage,
            final List<HlaSequenceLoci> winningSequences, final List<HlaSequenceLoci> winningNucSequences)
    {
        if(mConfig.TumorBam.isEmpty())
        {
            mTumorCoverage = HlaComplexCoverage.create(Lists.newArrayList());
            return;
        }

        List<Fragment> tumorNucleotideFrags = mNucleotideGeneEnrichment.enrich(mTumorBamReader.readFromBam());

        List<Fragment> tumorFragments = mAminoAcidPipeline.calcComparisonCoverageFragments(tumorNucleotideFrags);

        LL_LOGGER.info("calculating tumor coverage from frags({} highQual={})", tumorNucleotideFrags.size(), tumorFragments.size());

        List<FragmentAlleles> tumorFragAlleles = fragAlleleMapper.createFragmentAlleles(tumorFragments, winningSequences, winningNucSequences);

        fragAlleleMapper.checkHlaYSupport(mRefData.HlaYNucleotideSequences, tumorFragAlleles, tumorFragments, false, "tumor");

        mTumorCoverage = HlaComplexBuilder.calcProteinCoverage(tumorFragAlleles, winningAlleles);

        mTumorCoverage.populateMissingCoverage(winningAlleles);

        if(!mConfig.CopyNumberFile.isEmpty())
        {
            CopyNumberAssignment copyNumberAssignment = new CopyNumberAssignment();
            copyNumberAssignment.loadCopyNumberData(mConfig);

            copyNumberAssignment.assign(mConfig.Sample, winningAlleles, winningRefCoverage, mTumorCoverage, mTumorCopyNumber);
        }

        // SOMATIC VARIANTS
        SomaticVariantAnnotation variantAnnotation = new SomaticVariantAnnotation(
                mConfig, mRefData.HlaTranscriptData, mRefData.LociPositionFinder);

        if(variantAnnotation.getSomaticVariants().isEmpty())
            return;

        mSomaticVariants.addAll(variantAnnotation.getSomaticVariants());

        LL_LOGGER.info("calculating somatic variant allele coverage");

        boolean hasVcfData = mSomaticVariants.stream().anyMatch(x -> x.Context != null);
        LilacVCF lilacVCF = null;

        if(hasVcfData)
        {
            String vcfFilename = mConfig.outputPrefix() + ".lilac.somatic.vcf.gz";
            final VersionInfo version = new VersionInfo("lilac.version");
            lilacVCF = new LilacVCF(vcfFilename, mConfig.SomaticVariantsFile).writeHeader(version.toString());
        }

        for(SomaticVariant variant : mSomaticVariants)
        {
            List<HlaAlleleCoverage> variantCoverage = variantAnnotation.assignAlleleCoverage(variant, mTumorBamReader, winningSequences);
            List<HlaAllele> variantAlleles = variantCoverage.stream().map(x -> x.Allele).collect(Collectors.toList());
            LL_LOGGER.info("  {} -> {}}", variant, variantCoverage);

            if(lilacVCF != null && variant.Context != null)
                lilacVCF.writeVariant(variant.Context, variantAlleles);

            addVariant(mSomaticCodingCounts, variant, variantAlleles);
        }

        if(lilacVCF != null)
            lilacVCF.close();
    }

    public void extractRnaCoverage(
            final FragmentAlleleMapper fragAlleleMapper, final List<HlaAllele> winningAlleles,
            final List<HlaSequenceLoci> winningSequences, final List<HlaSequenceLoci> winningNucSequences)
    {
        if(mConfig.RnaBam.isEmpty())
        {
            mRnaCoverage = HlaComplexCoverage.create(Lists.newArrayList());
            return;
        }

        BamRecordReader rnaBamReader = new BamRecordReader(mConfig.RnaBam, mConfig.RefGenome, mRefData.HlaTranscriptData, mNucleotideFragFactory);

        List<Fragment> rnaNucleotideFrags = mNucleotideGeneEnrichment.enrich(rnaBamReader.readFromBam());

        List<Fragment> rnaFragments = mAminoAcidPipeline.calcComparisonCoverageFragments(rnaNucleotideFrags);

        LL_LOGGER.info("calculating RNA coverage from frags({} highQual={})", rnaNucleotideFrags.size(), rnaFragments.size());

        List<FragmentAlleles> rnaFragAlleles = fragAlleleMapper.createFragmentAlleles(rnaFragments, winningSequences, winningNucSequences);

        fragAlleleMapper.checkHlaYSupport(mRefData.HlaYNucleotideSequences, rnaFragAlleles, rnaFragments, false, "rna");

        mRnaCoverage = HlaComplexBuilder.calcProteinCoverage(rnaFragAlleles, winningAlleles);

        mRnaCoverage.populateMissingCoverage(winningAlleles);
    }

    public void writeFileOutputs()
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        mSummaryMetrics.log(mConfig.Sample);

        LL_LOGGER.info("writing output to {}", mConfig.OutputDir);
        String outputFile = mConfig.outputPrefix() + ".lilac.csv";
        String outputQCFile = mConfig.outputPrefix() + ".lilac.qc.csv";

        mSolutionSummary.write(outputFile);
        mSummaryMetrics.writefile(outputQCFile);

        HlaComplexFile.writeToFile(String.format("%s.candidates.coverage.csv", mConfig.outputPrefix()), mRankedComplexes);
        HlaComplexFile.writeFragmentAssignment(
                String.format("%s.candidates.fragments.csv", mConfig.outputPrefix()), mRankedComplexes, mRefFragAlleles);

        mAminoAcidPipeline.writeCounts(mConfig.outputPrefix());

        mRefAminoAcidCounts.writeVertically(String.format("%s.candidates.aminoacids.txt", mConfig.outputPrefix()));
        mRefNucleotideCounts.writeVertically(String.format("%s.candidates.nucleotides.txt", mConfig.outputPrefix()));

        FragmentUtils.writeFragmentData(String.format("%s.fragments.csv", mConfig.outputPrefix()), mRefNucleotideFrags);
    }

    private boolean validateFragments(final List<Fragment> fragments)
    {
        if(!mConfig.RunValidation)
            return true;

        List<Fragment> invalidFragments = fragments.stream().filter(x -> !x.validate()).collect(Collectors.toList());
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
            return  true;

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

    public void writeDatabaseResults(final CommandLine cmd)
    {
        if(!DatabaseAccess.hasDatabaseConfig(cmd))
            return;

        /*
        LL_LOGGER.info("Writing output to DB");
        DatabaseAccess dbAccess = DatabaseAccess.databaseAccess((CommandLine) cmd, (boolean) true);
        HlaType type = HlaFiles.type((String) outputFile, (String) outputQCFile);
        List typeDetails = HlaFiles.typeDetails((String) outputFile);
        dbAccess.writeHla(sample, type, typeDetails);
        */
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = LilacConfig.createOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }
        else if(cmd.hasOption(LOG_LEVEL))
        {
            Configurator.setRootLevel(Level.valueOf(cmd.getOptionValue(LOG_LEVEL)));
        }

        LilacApplication lilac = new LilacApplication(new LilacConfig(cmd));

        long startTime = System.currentTimeMillis();

        lilac.run();
        lilac.writeFileOutputs();

        if(DatabaseAccess.hasDatabaseConfig(cmd))
            lilac.writeDatabaseResults(cmd);

        long endTime = System.currentTimeMillis();
        double runTime = (endTime - startTime) / 1000.0;

        LL_LOGGER.info("Lilac complete, run time({}s)", String.format("%.2f", runTime));
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
