package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConfig.LOG_DEBUG;
import static com.hartwig.hmftools.lilac.LilacConfig.LOG_LEVEL;
import static com.hartwig.hmftools.lilac.LilacConstants.ALL_NUCLEOTIDE_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.LilacConstants.A_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.LilacConstants.B_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.LilacConstants.C_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_A;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_B;
import static com.hartwig.hmftools.lilac.LilacConstants.GENE_C;
import static com.hartwig.hmftools.lilac.coverage.CoverageCalcTask.proteinCoverage;
import static com.hartwig.hmftools.lilac.fragment.AminoAcidFragment.nucFragments;
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
import com.hartwig.hmftools.lilac.evidence.PhasedEvidenceValidation;
import com.hartwig.hmftools.lilac.fragment.AminoAcidFragment;
import com.hartwig.hmftools.lilac.fragment.AminoAcidFragmentPipeline;
import com.hartwig.hmftools.lilac.fragment.NucleotideGeneEnrichment;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaContext;
import com.hartwig.hmftools.lilac.hla.HlaContextFactory;
import com.hartwig.hmftools.lilac.fragment.NucleotideFragment;
import com.hartwig.hmftools.lilac.fragment.NucleotideFragmentFactory;
import com.hartwig.hmftools.lilac.hla.HlaCopyNumber;
import com.hartwig.hmftools.lilac.qc.AminoAcidQC;
import com.hartwig.hmftools.lilac.qc.BamQC;
import com.hartwig.hmftools.lilac.qc.CoverageQC;
import com.hartwig.hmftools.lilac.qc.HaplotypeQC;
import com.hartwig.hmftools.lilac.qc.HlaOut;
import com.hartwig.hmftools.lilac.qc.LilacQC;
import com.hartwig.hmftools.lilac.qc.SomaticVariantQC;
import com.hartwig.hmftools.lilac.read.FragmentAlleles;
import com.hartwig.hmftools.lilac.read.SAMRecordReader;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci;
import com.hartwig.hmftools.lilac.seq.HlaSequenceLociFile;
import com.hartwig.hmftools.lilac.variant.LilacVCF;
import com.hartwig.hmftools.lilac.variant.SomaticAlleleCoverage;
import com.hartwig.hmftools.lilac.variant.SomaticCodingCount;
import com.hartwig.hmftools.lilac.variant.SomaticVariants;

import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;
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
        LL_LOGGER.info("Starting LILAC with parameters:");
        mConfig.logParams();

        HlaContextFactory
                hlaContextFactory = new HlaContextFactory(A_EXON_BOUNDARIES, B_EXON_BOUNDARIES, C_EXON_BOUNDARIES);
        HlaContext hlaAContext = hlaContextFactory.hlaA();
        HlaContext hlaBContext = hlaContextFactory.hlaB();
        HlaContext hlaCContext = hlaContextFactory.hlaC();

        if(!mRefData.load(true))
        {
            LL_LOGGER.error("reference data loading failed");
            System.exit(1);
        }

        LL_LOGGER.info("Querying records from reference bam " + mConfig.ReferenceBam);

        NucleotideFragmentFactory nucleotideFragmentFactory = new NucleotideFragmentFactory(
                mConfig.MinBaseQual, mRefData.AminoAcidSequencesWithInserts, mRefData.AminoAcidSequencesWithDeletes,
                mRefData.LociPositionFinder);

        SAMRecordReader tumorBamReader =
                new SAMRecordReader(mConfig.TumorBam, mConfig.RefGenome, mRefData.HlaTranscriptData, nucleotideFragmentFactory);

        SAMRecordReader referenceBamReader =
                new SAMRecordReader(mConfig.ReferenceBam, mConfig.RefGenome, mRefData.HlaTranscriptData, nucleotideFragmentFactory);

        final List<NucleotideFragment> referenceNucleotideFragments = mNucleotideGeneEnrichment.enrich(referenceBamReader.readFromBam());
        final List<NucleotideFragment> tumorNucleotideFragments = Lists.newArrayList();
        if(!mConfig.TumorBam.isEmpty())
        {
            LL_LOGGER.info("Querying records from tumor bam " + mConfig.TumorBam);
            tumorNucleotideFragments.addAll(mNucleotideGeneEnrichment.enrich(tumorBamReader.readFromBam()));
        }

        validateNucleotideFragments(referenceNucleotideFragments);
        validateNucleotideFragments(tumorNucleotideFragments);

        AminoAcidFragmentPipeline aminoAcidPipeline = new AminoAcidFragmentPipeline(
                mConfig, referenceNucleotideFragments, tumorNucleotideFragments);

        // Enrich bam records
        List<AminoAcidFragment> aCandidateFragments = aminoAcidPipeline.referencePhasingFragments(hlaAContext);
        List<AminoAcidFragment> bCandidateFragments = aminoAcidPipeline.referencePhasingFragments(hlaBContext);
        List<AminoAcidFragment> cCandidateFragments = aminoAcidPipeline.referencePhasingFragments(hlaCContext);

        // Un-phased Candidates
        Candidates candidateFactory = new Candidates(mConfig, mRefData.NucleotideSequences, mRefData.AminoAcidSequences);
        List<HlaAllele> aUnphasedCandidates = candidateFactory.unphasedCandidates(hlaAContext, aCandidateFragments);
        List<HlaAllele> bUnphasedCandidates = candidateFactory.unphasedCandidates(hlaBContext, bCandidateFragments);
        List<HlaAllele> cUnphasedCandidates = candidateFactory.unphasedCandidates(hlaCContext, cCandidateFragments);
        List<HlaAllele> allUnphasedCandidates = Lists.newArrayList();
        allUnphasedCandidates.addAll(aUnphasedCandidates);
        allUnphasedCandidates.addAll(bUnphasedCandidates);
        allUnphasedCandidates.addAll(cUnphasedCandidates);

        // Phasing
        PhasedEvidenceFactory phasedEvidenceFactory = new PhasedEvidenceFactory(mConfig);
        List<PhasedEvidence> aPhasedEvidence = phasedEvidenceFactory.evidence(hlaAContext, aCandidateFragments);
        List<PhasedEvidence> bPhasedEvidence = phasedEvidenceFactory.evidence(hlaBContext, bCandidateFragments);
        List<PhasedEvidence> cPhasedEvidence = phasedEvidenceFactory.evidence(hlaCContext, cCandidateFragments);

        // validate phasing against expected sequences
        List<HlaSequenceLoci> expectedSequences = Lists.newArrayList();
        if(!mConfig.ExpectedAlleles.isEmpty())
        {
            expectedSequences.addAll(mRefData.AminoAcidSequences.stream()
                    .filter(x -> mConfig.ExpectedAlleles.contains(x.getAllele().asFourDigit())).collect(Collectors.toList()));

            PhasedEvidenceValidation.logInconsistentEvidence(GENE_A, aPhasedEvidence, expectedSequences);
            PhasedEvidenceValidation.logInconsistentEvidence(GENE_B, bPhasedEvidence, expectedSequences);
            PhasedEvidenceValidation.logInconsistentEvidence(GENE_C, cPhasedEvidence, expectedSequences);
        }

        // gather up all phased candidates
        List<HlaAllele> phasedCandidateAlleles = Lists.newArrayList();
        List<HlaAllele> aCandidates = candidateFactory.phasedCandidates(hlaAContext, aUnphasedCandidates, aPhasedEvidence);
        List<HlaAllele> bCandidates = candidateFactory.phasedCandidates(hlaBContext, bUnphasedCandidates, bPhasedEvidence);
        List<HlaAllele> cCandidates = candidateFactory.phasedCandidates(hlaCContext, cUnphasedCandidates, cPhasedEvidence);
        phasedCandidateAlleles.addAll(aCandidates);
        phasedCandidateAlleles.addAll(bCandidates);
        phasedCandidateAlleles.addAll(cCandidates);

        List<HlaAllele> candidateAlleles = phasedCandidateAlleles.stream().collect(Collectors.toList());
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
                    // .filter(x -> contains(allUnphasedCandidates, x)) // was previously only added if supported but now always
                    .collect(Collectors.toList());

            if(!missedCommonAlleles.isEmpty())
            {
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
                .filter(x -> candidateAlleles.contains(x.getAllele())).collect(Collectors.toList());

        // Coverage
        List<AminoAcidFragment> referenceCoverageFragments = aminoAcidPipeline.referenceCoverageFragments();
        validateAminoAcidFragments(referenceCoverageFragments);

        SequenceCount referenceAminoAcidCounts = SequenceCount.aminoAcids(mConfig.MinEvidence, referenceCoverageFragments);
        List<Integer> referenceAminoAcidHeterozygousLoci = referenceAminoAcidCounts.heterozygousLoci();
        SequenceCount referenceNucleotideCounts = SequenceCount.nucleotides(mConfig.MinEvidence, nucFragments(referenceCoverageFragments));
        List<Integer> referenceNucleotideHeterozygousLoci = referenceNucleotideCounts.heterozygousLoci().stream()
                .filter(x -> ALL_NUCLEOTIDE_EXON_BOUNDARIES.contains(x)).collect(Collectors.toList());

        List<HlaAllele> candidateAlleleSpecificProteins = candidateAlleles.stream().map(x -> x.asFourDigit()).collect(Collectors.toList());
        List<HlaSequenceLoci> candidateAminoAcidSequences = mRefData.AminoAcidSequences.stream()
                .filter(x -> candidateAlleles.contains(x.getAllele())).collect(Collectors.toList());

        List<HlaSequenceLoci> candidateNucleotideSequences = mRefData.NucleotideSequences.stream()
            .filter(x -> candidateAlleleSpecificProteins.contains(x.getAllele().asFourDigit())).collect(Collectors.toList());

        List<FragmentAlleles> referenceFragmentAlleles = FragmentAlleles.create(
                referenceCoverageFragments, referenceAminoAcidHeterozygousLoci, candidateAminoAcidSequences,
                referenceNucleotideHeterozygousLoci, candidateNucleotideSequences);

        FragmentAlleles.applyUniqueStopLossFragments(
                referenceFragmentAlleles, referenceBamReader.stopLossOnCIndels(), mRefData.StopLossRecoveryAlleles);

        // build and score complexes
        HlaComplexBuilder complexBuilder = new HlaComplexBuilder(mConfig, mRefData);
        List<HlaComplex> complexes = complexBuilder.buildComplexes(referenceFragmentAlleles, candidateAlleles);

        LL_LOGGER.info("calculating coverage of {} complexes", complexes.size());
        ComplexCoverageCalculator complexCalculator = new ComplexCoverageCalculator(mConfig);
        List<HlaComplexCoverage> calculatedComplexes = complexCalculator.calculateComplexCoverages(referenceFragmentAlleles, complexes);

        HlaComplexCoverageRanking complexRanker = new HlaComplexCoverageRanking(mConfig.MaxDistanceFromTopScore, mRefData);
        List<HlaComplexCoverage> referenceRankedComplexes = complexRanker.rankCandidates(calculatedComplexes);

        if(referenceRankedComplexes.isEmpty())
        {
            LL_LOGGER.fatal("failed to calculate complex coverage");
            System.exit(1);
        }

        if(!expectedSequences.isEmpty())
        {
            HlaComplexCoverage expectedCoverage = proteinCoverage(
                    referenceFragmentAlleles, expectedSequences.stream().map(x -> x.getAllele()).collect(Collectors.toList()));
            LL_LOGGER.info("expected allele coverage: {}", expectedCoverage);
        }

        HlaComplexCoverage winningReferenceCoverage = referenceRankedComplexes.get(0).expandToSixAlleles();
        List<HlaAllele> winningAlleles = winningReferenceCoverage.getAlleleCoverage().stream().map(x -> x.Allele).collect(Collectors.toList());
        List<HlaSequenceLoci> winningSequences = candidateSequences.stream()
                .filter(x -> winningAlleles.contains(x.getAllele())).collect(Collectors.toList());

        LL_LOGGER.info("{}", HlaComplexFile.header());

        for (HlaComplexCoverage rankedComplex : referenceRankedComplexes)
        {
            LL_LOGGER.info(HlaComplexFile.asString(rankedComplex));
        }

        LL_LOGGER.info("{} - REF - {} CANDIDATES, WINNING ALLELES: {}",
                mConfig.Sample, referenceRankedComplexes.size(),
                HlaAllele.toString(winningReferenceCoverage.getAlleleCoverage().stream().map(x -> x.Allele).collect(Collectors.toList())));

        HlaComplexCoverage winningTumorCoverage = null;
        List<HlaCopyNumber> winningTumorCopyNumber = null;
        List<VariantContextDecorator> somaticVariants = Lists.newArrayList();
        List<SomaticCodingCount> somaticCodingCounts = SomaticCodingCount.create(winningAlleles);

        if(!mConfig.TumorBam.isEmpty())
        {
            LL_LOGGER.info("Calculating tumor coverage of winning alleles");

            List<FragmentAlleles> tumorFragmentAlleles = FragmentAlleles.create(
                    aminoAcidPipeline.tumorCoverageFragments(),
                    referenceAminoAcidHeterozygousLoci, candidateAminoAcidSequences,
                    referenceNucleotideHeterozygousLoci, candidateNucleotideSequences);

            winningTumorCoverage = proteinCoverage(tumorFragmentAlleles, winningAlleles).expandToSixAlleles();

            LL_LOGGER.info("Calculating tumor copy number of winning alleles");
            winningTumorCopyNumber = HlaCopyNumber.alleleCopyNumber(winningAlleles, mConfig.GeneCopyNumberFile, winningTumorCoverage);

            // SOMATIC VARIANTS
            somaticVariants.addAll(new SomaticVariants(mConfig, mRefData.HlaTranscriptData).readSomaticVariants());

            if(!somaticVariants.isEmpty())
            {
                LL_LOGGER.info("Calculating somatic variant allele coverage");

                String vcfFilename = mConfig.outputPrefix() + ".lilac.somatic.vcf.gz";
                LilacVCF lilacVCF = new LilacVCF(vcfFilename, mConfig.SomaticVcf).writeHeader(version.toString());

                SomaticAlleleCoverage somaticCoverageFactory = new SomaticAlleleCoverage(
                        mConfig, referenceAminoAcidHeterozygousLoci, mRefData.LociPositionFinder, somaticVariants, winningSequences);

                for (VariantContextDecorator variant : somaticVariants)
                {
                    List<HlaAlleleCoverage> variantCoverage = somaticCoverageFactory.alleleCoverage(variant, tumorBamReader);
                    List<HlaAllele> variantAlleles = variantCoverage.stream().map(x -> x.Allele).collect(Collectors.toList());
                    LL_LOGGER.info("    {} -> {}}", variant, variantCoverage);
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

        HlaOut output = HlaOut.create(winningReferenceCoverage, winningTumorCoverage, winningTumorCopyNumber, somaticCodingCounts);

        // QC
        LL_LOGGER.info("Calculating QC Statistics");
        SomaticVariantQC somaticVariantQC = SomaticVariantQC.create(somaticVariants.size(), somaticCodingCounts);
        AminoAcidQC aminoAcidQC = AminoAcidQC.create(winningSequences, referenceAminoAcidCounts);

        List<PhasedEvidence> combinedPhasedEvidence = Lists.newArrayList();
        combinedPhasedEvidence.addAll(aPhasedEvidence);
        combinedPhasedEvidence.addAll(bPhasedEvidence);
        combinedPhasedEvidence.addAll(cPhasedEvidence);
        HaplotypeQC haplotypeQC = HaplotypeQC.create(3, winningSequences, combinedPhasedEvidence, referenceAminoAcidCounts);

        BamQC bamQC = BamQC.create(referenceBamReader);
        CoverageQC coverageQC = CoverageQC.create(referenceNucleotideFragments.size(), winningReferenceCoverage);
        LilacQC lilacQC = LilacQC.create(aminoAcidQC, bamQC, coverageQC, haplotypeQC, somaticVariantQC);

        LL_LOGGER.info("QC Stats:");
        LL_LOGGER.info("  {}", lilacQC.header());
        LL_LOGGER.info("  {}", lilacQC.body());

        LL_LOGGER.info("Writing output to {}", mConfig.OutputDir);
        String outputFile = mConfig.outputPrefix() + ".lilac.txt";
        String outputQCFile = mConfig.outputPrefix() + ".lilac.qc.txt";

        output.write(outputFile);
        lilacQC.writefile(outputQCFile);

        List<HlaSequenceLoci> candidatesToWrite = candidateSequences.stream().collect(Collectors.toList());
        candidatesToWrite.addAll(expectedSequences);

        if(mRefData.getDeflatedSequenceTemplate() != null)
            candidatesToWrite.add(mRefData.getDeflatedSequenceTemplate());

        HlaSequenceLociFile.write(
                String.format("%s.candidates.sequences.txt", mConfig.outputPrefix()),
                A_EXON_BOUNDARIES, B_EXON_BOUNDARIES, C_EXON_BOUNDARIES, candidatesToWrite);

        HlaComplexFile.writeToFile(referenceRankedComplexes,
                String.format("%s.candidates.coverage.txt", mConfig.outputPrefix()));

        referenceAminoAcidCounts.writeVertically(String.format("%s.candidates.aminoacids.txt", mConfig.outputPrefix()));
        referenceNucleotideCounts.writeVertically(String.format("%s.candidates.nucleotides.txt", mConfig.outputPrefix()));

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
    }

    private boolean validateAminoAcidFragments(final List<AminoAcidFragment> fragments)
    {
        if(!mConfig.RunValidation)
            return true;

        List<AminoAcidFragment> invalidFragments = fragments.stream().filter(x -> !x.validate()).collect(Collectors.toList());
        if(invalidFragments.isEmpty())
            return true;


        LL_LOGGER.warn("has {} invalid amino-acid fragments", invalidFragments.size());
        return false;
    }

    private boolean validateNucleotideFragments(final List<NucleotideFragment> fragments)
    {
        if(!mConfig.RunValidation)
            return true;

        List<NucleotideFragment> invalidFragments = fragments.stream().filter(x -> !x.validate()).collect(Collectors.toList());
        if(invalidFragments.isEmpty())
            return true;


        LL_LOGGER.warn("has {} invalid nucleotide fragments", invalidFragments.size());
        return false;
    }

    @Override
    public void close()
    {
        LL_LOGGER.info("Finished in " + (System.currentTimeMillis() - mStartTime) / (long) 1000 + " seconds");
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
