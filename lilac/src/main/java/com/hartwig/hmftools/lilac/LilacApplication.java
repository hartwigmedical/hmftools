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
import static com.hartwig.hmftools.lilac.LilacConstants.LOG_UNMATCHED_HAPLOTYPE_SUPPORT;
import static com.hartwig.hmftools.lilac.candidates.NucleotideFiltering.calcNucleotideHeterogygousLoci;
import static com.hartwig.hmftools.lilac.coverage.FragmentAlleles.createFragmentAlleles;
import static com.hartwig.hmftools.lilac.variant.SomaticCodingCount.addVariant;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.lilac.candidates.Candidates;
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
import com.hartwig.hmftools.lilac.fragment.NucleotideGeneEnrichment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaContext;
import com.hartwig.hmftools.lilac.hla.HlaContextFactory;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.fragment.NucleotideFragmentFactory;
import com.hartwig.hmftools.lilac.variant.HlaCopyNumber;
import com.hartwig.hmftools.lilac.qc.AminoAcidQC;
import com.hartwig.hmftools.lilac.qc.BamQC;
import com.hartwig.hmftools.lilac.qc.CoverageQC;
import com.hartwig.hmftools.lilac.qc.HaplotypeQC;
import com.hartwig.hmftools.lilac.qc.HlaOut;
import com.hartwig.hmftools.lilac.qc.LilacQC;
import com.hartwig.hmftools.lilac.qc.SomaticVariantQC;
import com.hartwig.hmftools.lilac.coverage.FragmentAlleles;
import com.hartwig.hmftools.lilac.read.SAMRecordReader;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.variant.LilacVCF;
import com.hartwig.hmftools.lilac.variant.SomaticAlleleCoverage;
import com.hartwig.hmftools.lilac.variant.SomaticCodingCount;
import com.hartwig.hmftools.lilac.variant.SomaticVariantFinder;

import java.util.Collections;
import java.util.List;
import java.util.Map;
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

public class LilacApplication implements AutoCloseable, Runnable
{
    private final long mStartTime;
    private final NucleotideGeneEnrichment mNucleotideGeneEnrichment;
    private final LilacConfig mConfig;
    private final ReferenceData mRefData;

    public LilacApplication(final CommandLine cmd, final LilacConfig config)
    {
        mConfig = config;
        mRefData = new ReferenceData(mConfig.ResourceDir, mConfig);

        mStartTime = System.currentTimeMillis();

        mNucleotideGeneEnrichment = new NucleotideGeneEnrichment(A_EXON_BOUNDARIES, B_EXON_BOUNDARIES, C_EXON_BOUNDARIES);
    }

    public void run()
    {
        if(!mConfig.isValid())
        {
            LL_LOGGER.warn("invalid config, exiting");
            System.exit(1);
        }

        VersionInfo version = new VersionInfo("lilac.version");
        LL_LOGGER.info("starting Lilac with parameters:");
        mConfig.logParams();

        HlaContextFactory
                hlaContextFactory = new HlaContextFactory(A_EXON_BOUNDARIES, B_EXON_BOUNDARIES, C_EXON_BOUNDARIES);
        HlaContext hlaAContext = hlaContextFactory.hlaA();
        HlaContext hlaBContext = hlaContextFactory.hlaB();
        HlaContext hlaCContext = hlaContextFactory.hlaC();

        if(!mRefData.load())
        {
            LL_LOGGER.error("reference data loading failed");
            System.exit(1);
        }

        boolean allValid = true;

        LL_LOGGER.info("querying records from reference bam " + mConfig.ReferenceBam);

        NucleotideFragmentFactory nucleotideFragmentFactory = new NucleotideFragmentFactory(
                mConfig.MinBaseQual, mRefData.AminoAcidSequencesWithInserts, mRefData.AminoAcidSequencesWithDeletes,
                mRefData.LociPositionFinder);

        SAMRecordReader tumorBamReader =
                new SAMRecordReader(mConfig.TumorBam, mConfig.RefGenome, mRefData.HlaTranscriptData, nucleotideFragmentFactory);

        SAMRecordReader referenceBamReader =
                new SAMRecordReader(mConfig.ReferenceBam, mConfig.RefGenome, mRefData.HlaTranscriptData, nucleotideFragmentFactory);

        final List<Fragment> refNucleotideFrags = mNucleotideGeneEnrichment.enrich(referenceBamReader.readFromBam());
        final List<Fragment> tumorNucleotideFrags = Lists.newArrayList();
        if(!mConfig.TumorBam.isEmpty())
        {
            LL_LOGGER.info("querying records from tumor bam " + mConfig.TumorBam);
            tumorNucleotideFrags.addAll(mNucleotideGeneEnrichment.enrich(tumorBamReader.readFromBam()));
        }

        allValid &= validateFragments(refNucleotideFrags);
        allValid &= validateFragments(tumorNucleotideFrags);

        AminoAcidFragmentPipeline aminoAcidPipeline = new AminoAcidFragmentPipeline(mConfig, refNucleotideFrags, tumorNucleotideFrags);

        /*
        // repeat the heterozygous-loci check on all candidates
        List<Fragment> refAminoAcidFragments = aminoAcidPipeline.referenceAminoAcidFragments();
        SequenceCount refAminoAcidCounts = SequenceCount.aminoAcids(mConfig.MinEvidence, refAminoAcidFragments);
        List<HlaSequenceLoci> supportedAminoAcidSequences = filterCandidatesByAminoAcidLoci(mRefData.AminoAcidSequences, refAminoAcidCounts);

        List<HlaSequenceLoci> supportedAminoAcidSequences = mRefData.AminoAcidSequences;

        LL_LOGGER.info("filtered amino acid sequences {} -> {} by min support",
                mRefData.AminoAcidSequences.size(), supportedAminoAcidSequences.size());
        */

        // apply special filtering and splice checks on fragments, just for use in phasing
        List<Fragment> aCandidateFrags = aminoAcidPipeline.referencePhasingFragments(hlaAContext);
        List<Fragment> bCandidateFrags = aminoAcidPipeline.referencePhasingFragments(hlaBContext);
        List<Fragment> cCandidateFrags = aminoAcidPipeline.referencePhasingFragments(hlaCContext);

        // determine un-phased Candidates
        Candidates candidateFactory = new Candidates(mConfig, mRefData.NucleotideSequences, mRefData.AminoAcidSequences);
        List<HlaAllele> aUnphasedCandidates = candidateFactory.unphasedCandidates(hlaAContext, aCandidateFrags);
        List<HlaAllele> bUnphasedCandidates = candidateFactory.unphasedCandidates(hlaBContext, bCandidateFrags);
        List<HlaAllele> cUnphasedCandidates = candidateFactory.unphasedCandidates(hlaCContext, cCandidateFrags);

        // determine phasing of amino acids
        PhasedEvidenceFactory phasedEvidenceFactory = new PhasedEvidenceFactory(mConfig);
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
        candidateAlleles.addAll(aCandidates);
        candidateAlleles.addAll(bCandidates);
        candidateAlleles.addAll(cCandidates);

        List<HlaAllele> recoveredAlleles = Lists.newArrayList();

        if(!mConfig.RestrictedAlleles.isEmpty())
        {
            // if using a set of restricted allels, make no concession for any missed or common alleles
            mConfig.ExpectedAlleles.stream().filter(x -> !candidateAlleles.contains(x)).forEach(x -> candidateAlleles.add(x));
        }
        else
        {
            List<HlaAllele> missedStopLossAlleles = mRefData.StopLossRecoveryAlleles.stream()
                    .filter(x -> !candidateAlleles.contains(x)).collect(Collectors.toList());

            if(referenceBamReader.stopLossOnCIndels() > 0 && !missedStopLossAlleles.isEmpty())
            {
                recoveredAlleles.addAll(missedStopLossAlleles);
                LL_LOGGER.info("recovering {} stop loss alleles with {} fragments",
                        HlaAllele.toString(missedStopLossAlleles), referenceBamReader.stopLossOnCIndels());
            }

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
        List<Fragment> refAminoAcidFragments = aminoAcidPipeline.getReferenceFragments();
        SequenceCount refAminoAcidCounts = SequenceCount.aminoAcids(mConfig.MinEvidence, refAminoAcidFragments);
        SequenceCount refNucleotideCounts = SequenceCount.nucleotides(mConfig.MinEvidence, refAminoAcidFragments);

        //List<Integer> refAminoAcidHetLoci = refAminoAcidCounts.heterozygousLoci();

        List<Integer> refAminoAcidHetLoci = refAminoAcidCounts.heterozygousLoci();
        Map<String,List<Integer>> refNucleotideHetLoci = calcNucleotideHeterogygousLoci(refNucleotideCounts.heterozygousLoci());

        List<HlaAllele> candidateAlleleSpecificProteins = candidateAlleles.stream().map(x -> x.asFourDigit()).collect(Collectors.toList());

        List<HlaSequenceLoci> candidateNucleotideSequences = mRefData.NucleotideSequences.stream()
            .filter(x -> candidateAlleleSpecificProteins.contains(x.Allele.asFourDigit())).collect(Collectors.toList());

        List<FragmentAlleles> refFragAlleles = FragmentAlleles.createFragmentAlleles(
                refAminoAcidFragments, refAminoAcidHetLoci, candidateSequences, aminoAcidPipeline.getReferenceAminoAcids(),
                refNucleotideHetLoci, candidateNucleotideSequences, aminoAcidPipeline.getReferenceNucleotides());

        FragmentAlleles.applyUniqueStopLossFragments(
                refFragAlleles, referenceBamReader.stopLossOnCIndels(), mRefData.StopLossRecoveryAlleles);

        FragmentAlleles.checkHlaYSupport(
                mConfig.Sample, mRefData.HlaYNucleotideSequences, refAminoAcidHetLoci, refFragAlleles, refAminoAcidFragments);

        // build and score complexes
        HlaComplexBuilder complexBuilder = new HlaComplexBuilder(mConfig, mRefData);
        List<HlaComplex> complexes = complexBuilder.buildComplexes(
                refFragAlleles, candidateAlleles, recoveredAlleles, mRefData.getWildcardAlleles());

        LL_LOGGER.info("calculating coverage of {} complexes", complexes.size());
        ComplexCoverageCalculator complexCalculator = new ComplexCoverageCalculator(mConfig.Threads);
        List<HlaComplexCoverage> calculatedComplexes = complexCalculator.calculateComplexCoverages(refFragAlleles, complexes);

        HlaComplexCoverageRanking complexRanker = new HlaComplexCoverageRanking(mConfig.TopScoreThreshold, mRefData);
        List<HlaComplexCoverage> referenceRankedComplexes = complexRanker.rankCandidates(calculatedComplexes, recoveredAlleles);

        if(referenceRankedComplexes.isEmpty())
        {
            LL_LOGGER.fatal("failed to calculate complex coverage");
            System.exit(1);
        }

        if(!expectedSequences.isEmpty())
        {
            HlaComplexCoverage expectedCoverage = HlaComplexBuilder.calcProteinCoverage(
                    refFragAlleles, expectedSequences.stream().map(x -> x.Allele).collect(Collectors.toList()));
            LL_LOGGER.info("expected allele coverage: {}", expectedCoverage);
        }

        HlaComplexCoverage winningRefCoverage = referenceRankedComplexes.get(0).expandToSixAlleles();
        List<HlaAllele> winningAlleles = winningRefCoverage.getAlleleCoverage().stream().map(x -> x.Allele).collect(Collectors.toList());
        List<HlaSequenceLoci> winningSequences = candidateSequences.stream()
                .filter(x -> winningAlleles.contains(x.Allele)).collect(Collectors.toList());

        LL_LOGGER.info("{}", HlaComplexFile.header());

        for (HlaComplexCoverage rankedComplex : referenceRankedComplexes)
        {
            LL_LOGGER.info(HlaComplexFile.asString(rankedComplex));
        }

        // log key results for fast post-run analysis
        StringJoiner totalCoverages = new StringJoiner(",");
        winningRefCoverage.getAlleleCoverage().forEach(x -> totalCoverages.add(String.format("%.0f",x.TotalCoverage)));

        LL_LOGGER.info("WINNERS: {},{},{},{}",
                mConfig.Sample, referenceRankedComplexes.size(), HlaAllele.toString(winningRefCoverage.getAlleles()), totalCoverages);

        HlaComplexCoverage winningTumorCoverage = null;
        List<HlaCopyNumber> winningTumorCopyNumber = null;
        List<VariantContextDecorator> somaticVariants = Lists.newArrayList();
        List<SomaticCodingCount> somaticCodingCounts = SomaticCodingCount.create(winningAlleles);

        if(!mConfig.TumorBam.isEmpty())
        {
            LL_LOGGER.info("calculating tumor coverage of winning alleles");

            List<FragmentAlleles> tumorFragAlleles = createFragmentAlleles(
                    aminoAcidPipeline.tumorCoverageFragments(), refAminoAcidHetLoci, candidateSequences,
                    aminoAcidPipeline.getReferenceAminoAcids(),
                    refNucleotideHetLoci, candidateNucleotideSequences, aminoAcidPipeline.getReferenceNucleotides());

            winningTumorCoverage = HlaComplexBuilder.calcProteinCoverage(tumorFragAlleles, winningAlleles).expandToSixAlleles();

            if(!mConfig.GeneCopyNumberFile.isEmpty())
            {
                LL_LOGGER.info("calculating tumor copy number of winning alleles");
                winningTumorCopyNumber = HlaCopyNumber.alleleCopyNumber(winningAlleles, mConfig.GeneCopyNumberFile, winningTumorCoverage);
            }

            // SOMATIC VARIANTS
            SomaticVariantFinder somaticVariantFinder = new SomaticVariantFinder(mConfig, mRefData.HlaTranscriptData);
            somaticVariants.addAll(somaticVariantFinder.readSomaticVariants());

            if(!somaticVariants.isEmpty())
            {
                LL_LOGGER.info("calculating somatic variant allele coverage");

                String vcfFilename = mConfig.outputPrefix() + ".lilac.somatic.vcf.gz";
                LilacVCF lilacVCF = new LilacVCF(vcfFilename, mConfig.SomaticVcf).writeHeader(version.toString());

                SomaticAlleleCoverage somaticCoverageFactory = new SomaticAlleleCoverage(
                        mConfig, refAminoAcidHetLoci, mRefData.LociPositionFinder, somaticVariants, winningSequences);

                for(VariantContextDecorator variant : somaticVariants)
                {
                    List<HlaAlleleCoverage> variantCoverage = somaticCoverageFactory.alleleCoverage(variant, tumorBamReader);
                    List<HlaAllele> variantAlleles = variantCoverage.stream().map(x -> x.Allele).collect(Collectors.toList());
                    LL_LOGGER.info("  {} -> {}}", variant, variantCoverage);
                    lilacVCF.writeVariant(variant.context(), variantAlleles);
                    addVariant(somaticCodingCounts, variant, variantAlleles);
                }
            }
        }
        else
        {
            winningTumorCoverage = HlaComplexCoverage.create(Lists.newArrayList());
            winningTumorCopyNumber = HlaCopyNumber.alleleCopyNumber(winningAlleles);
        }

        HlaOut output = HlaOut.create(winningRefCoverage, winningTumorCoverage, winningTumorCopyNumber, somaticCodingCounts);

        // create various QC and other metrics
        LL_LOGGER.info("calculating QC Statistics");
        SomaticVariantQC somaticVariantQC = SomaticVariantQC.create(somaticVariants.size(), somaticCodingCounts);

        List<PhasedEvidence> combinedPhasedEvidence = Lists.newArrayList();
        combinedPhasedEvidence.addAll(aPhasedEvidence);
        combinedPhasedEvidence.addAll(bPhasedEvidence);
        combinedPhasedEvidence.addAll(cPhasedEvidence);
        HaplotypeQC haplotypeQC = HaplotypeQC.create(
                LOG_UNMATCHED_HAPLOTYPE_SUPPORT, winningSequences, combinedPhasedEvidence, refAminoAcidCounts);

        AminoAcidQC aminoAcidQC = AminoAcidQC.create(winningSequences, refAminoAcidCounts, haplotypeQC.UnmatchedHaplotypes);

        BamQC bamQC = BamQC.create(referenceBamReader);
        CoverageQC coverageQC = CoverageQC.create(refNucleotideFrags.size(), winningRefCoverage);
        LilacQC lilacQC = LilacQC.create(aminoAcidQC, bamQC, coverageQC, haplotypeQC, somaticVariantQC);

        LL_LOGGER.info("QC Stats:");
        LL_LOGGER.info("  {}", lilacQC.header());
        LL_LOGGER.info("  {}", lilacQC.body());

        LL_LOGGER.info("writing output to {}", mConfig.OutputDir);
        String outputFile = mConfig.outputPrefix() + ".lilac.tsv";
        String outputQCFile = mConfig.outputPrefix() + ".lilac.qc.tsv";

        output.write(outputFile);
        lilacQC.writefile(outputQCFile);

        HlaComplexFile.writeToFile(String.format("%s.candidates.coverage.tsv", mConfig.outputPrefix()), referenceRankedComplexes);
        HlaComplexFile.writeFragmentAssignment(
                String.format("%s.candidates.fragments.csv", mConfig.outputPrefix()), referenceRankedComplexes, refFragAlleles);

        aminoAcidPipeline.writeCounts(mConfig.outputPrefix());

        refAminoAcidCounts.writeVertically(String.format("%s.candidates.aminoacids.txt", mConfig.outputPrefix()));
        refNucleotideCounts.writeVertically(String.format("%s.candidates.nucleotides.txt", mConfig.outputPrefix()));

        /*
        if(DatabaseAccess.hasDatabaseConfig((CommandLine) cmd))
        {
            LL_LOGGER.info("Writing output to DB");
            DatabaseAccess dbAccess = DatabaseAccess.databaseAccess((CommandLine) cmd, (boolean) true);
            HlaType type = HlaFiles.type((String) outputFile, (String) outputQCFile);
            List typeDetails = HlaFiles.typeDetails((String) outputFile);
            dbAccess.writeHla(sample, type, typeDetails);
        }
         */

        if(!allValid)
        {
            LL_LOGGER.error("failed validation");
        }
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

    @Override
    public void close()
    {
        LL_LOGGER.info("run time: {} seconds", (System.currentTimeMillis() - mStartTime) / (long) 1000);
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final VersionInfo version = new VersionInfo("lilac.version");
        LL_LOGGER.info("Lilac version: {}", version.version());

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

        LilacApplication lilac = new LilacApplication(cmd, new LilacConfig(cmd));
        lilac.run();
        lilac.close();

        LL_LOGGER.info("Lilac complete");

    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
