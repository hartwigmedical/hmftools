package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConfig.LOG_DEBUG;
import static com.hartwig.hmftools.lilac.LilacConfig.LOG_LEVEL;
import static com.hartwig.hmftools.lilac.LilacConstants.ALL_NUCLEOTIDE_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.LilacConstants.A_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.LilacConstants.B_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.LilacConstants.C_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.LilacConstants.DEFLATE_TEMPLATE;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_TRANSCRIPTS;
import static com.hartwig.hmftools.lilac.LilacConstants.LOCI_POSITION;
import static com.hartwig.hmftools.lilac.coverage.CoverageCalcTask.proteinCoverage;
import static com.hartwig.hmftools.lilac.fragment.AminoAcidFragment.nucFragments;
import static com.hartwig.hmftools.lilac.hla.HlaAllele.contains;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.lilac.candidates.Candidates;
import com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage;
import com.hartwig.hmftools.lilac.coverage.HlaComplex;
import com.hartwig.hmftools.lilac.coverage.HlaComplexCoverage;
import com.hartwig.hmftools.lilac.coverage.HlaComplexCoverageFactory;
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
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;
import java.util.stream.Collectors;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class LilacApplication implements AutoCloseable, Runnable
{
    private final long mStartTime;
    private final ExecutorService mExecutorService;
    private final NucleotideGeneEnrichment mNucleotideGeneEnrichment;
    private final LilacConfig mConfig;
    private final ReferenceData mRefData;

    public LilacApplication(final CommandLine cmd, final LilacConfig config)
    {
        mConfig = config;
        mRefData = new ReferenceData(mConfig.ResourceDir, mConfig);

        mStartTime = System.currentTimeMillis();

        final ThreadFactory namedThreadFactory = new ThreadFactoryBuilder().setNameFormat("Lilac-%d").build();

        mExecutorService = Executors.newFixedThreadPool(mConfig.Threads, namedThreadFactory);
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

        HlaSequenceLoci deflatedSequenceTemplate = mRefData.AminoAcidSequences.stream()
                .filter(x -> x.getAllele().matches(DEFLATE_TEMPLATE)).findFirst().orElse(null);

        LL_LOGGER.info("Querying records from reference bam " + mConfig.ReferenceBam);

        NucleotideFragmentFactory nucleotideFragmentFactory = new NucleotideFragmentFactory(
                mConfig.MinBaseQual, mRefData.AminoAcidSequencesWithInserts, mRefData.AminoAcidSequencesWithDeletes, LOCI_POSITION);

        SAMRecordReader tumorBamReader =
                new SAMRecordReader(mConfig.TumorBam, mConfig.RefGenome, HLA_TRANSCRIPTS, nucleotideFragmentFactory);

        SAMRecordReader referenceBamReader =
                new SAMRecordReader(mConfig.ReferenceBam, mConfig.RefGenome, HLA_TRANSCRIPTS, nucleotideFragmentFactory);

        final List<NucleotideFragment> referenceNucleotideFragments = mNucleotideGeneEnrichment.enrich(referenceBamReader.readFromBam());
        final List<NucleotideFragment> tumorNucleotideFragments = Lists.newArrayList();
        if(!mConfig.TumorBam.isEmpty())
        {
            LL_LOGGER.info("Querying records from tumor bam " + mConfig.TumorBam);
            tumorNucleotideFragments.addAll(mNucleotideGeneEnrichment.enrich(tumorBamReader.readFromBam()));
        }

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

        // Validate phasing against expected sequences
        List<HlaSequenceLoci> expectedSequences = mRefData.AminoAcidSequences.stream()
                .filter(x -> contains(mConfig.ExpectedAlleles, x.getAllele().asFourDigit())).collect(Collectors.toList());
        
        PhasedEvidenceValidation.logInconsistentEvidence("A", aPhasedEvidence, expectedSequences);
        PhasedEvidenceValidation.logInconsistentEvidence("B", bPhasedEvidence, expectedSequences);
        PhasedEvidenceValidation.logInconsistentEvidence("C", cPhasedEvidence, expectedSequences);

        // Phased Candidates
        List<HlaAllele> aCandidates = candidateFactory.phasedCandidates(hlaAContext, aUnphasedCandidates, aPhasedEvidence);
        List<HlaAllele> bCandidates = candidateFactory.phasedCandidates(hlaBContext, bUnphasedCandidates, bPhasedEvidence);
        List<HlaAllele> cCandidates = candidateFactory.phasedCandidates(hlaCContext, cUnphasedCandidates, cPhasedEvidence);
        List<HlaAllele> phasedCandidateAlleles = Lists.newArrayList();
        phasedCandidateAlleles.addAll(aCandidates);
        phasedCandidateAlleles.addAll(bCandidates);
        phasedCandidateAlleles.addAll(cCandidates);

        List<HlaAllele> missingStopLossAlleles = mRefData.StopLossRecoveryAlleles.stream()
                .filter(x -> !phasedCandidateAlleles.contains(x)).collect(Collectors.toList());

        List<HlaAllele> stopLossRecovery = Lists.newArrayList();
        
        if(referenceBamReader.stopLossOnCIndels() > 0 && !missingStopLossAlleles.isEmpty()) 
        {
            stopLossRecovery.addAll(missingStopLossAlleles); 
            LL_LOGGER.info("Identified {} stop loss fragments", referenceBamReader.stopLossOnCIndels());
            LL_LOGGER.info("  recovered stop loss candidates: {}", HlaAllele.toString(missingStopLossAlleles));
        }

        // Common Candidate Recovery
        LL_LOGGER.info("Recovering common un-phased candidates:");

        List<HlaAllele> recoveredAlleles = mRefData.CommonAlleles.stream()
                .filter(x -> !contains(phasedCandidateAlleles, x))
                .filter(x -> contains(allUnphasedCandidates, x))
                .collect(Collectors.toList());

        recoveredAlleles.addAll(stopLossRecovery);

        List<HlaAllele> candidateAlleles = phasedCandidateAlleles.stream().collect(Collectors.toList());
        candidateAlleles.addAll(recoveredAlleles);

        List<HlaAllele> missingExpected = mConfig.ExpectedAlleles.stream().filter(x -> !candidateAlleles.contains(x)).collect(Collectors.toList());

        if (!missingExpected.isEmpty())
        {
            LL_LOGGER.info("  including {} expected alleles as candidates: {}",
                    missingExpected.size(), HlaAllele.toString(missingExpected));

            candidateAlleles.addAll(missingExpected);
        }

        List<HlaSequenceLoci> candidateSequences = mRefData.AminoAcidSequences.stream()
                .filter(x -> contains(candidateAlleles, x.getAllele())).collect(Collectors.toList());

        if (!recoveredAlleles.isEmpty())
        {
            LL_LOGGER.info("  recovered {} common candidate alleles: {}",
                    recoveredAlleles.size(), HlaAllele.toString(recoveredAlleles));
        }
        else
        {
            LL_LOGGER.info("  recovered 0 common candidate alleles");
        }

        // Coverage
        List<AminoAcidFragment> referenceCoverageFragments = aminoAcidPipeline.referenceCoverageFragments();
        SequenceCount referenceAminoAcidCounts = SequenceCount.aminoAcids(mConfig.MinEvidence, referenceCoverageFragments);
        List<Integer> referenceAminoAcidHeterozygousLoci = referenceAminoAcidCounts.heterozygousLoci();
        SequenceCount referenceNucleotideCounts = SequenceCount.nucleotides(mConfig.MinEvidence, nucFragments(referenceCoverageFragments));
        List<Integer> referenceNucleotideHeterozygousLoci = referenceNucleotideCounts.heterozygousLoci().stream()
                .filter(x -> ALL_NUCLEOTIDE_EXON_BOUNDARIES.contains(x)).collect(Collectors.toList());

        List<HlaAllele> candidateAlleleSpecificProteins = candidateAlleles.stream().map(x -> x.asFourDigit()).collect(Collectors.toList());
        List<HlaSequenceLoci> candidateAminoAcidSequences = mRefData.AminoAcidSequences.stream()
                .filter(x -> contains(candidateAlleles, x.getAllele())).collect(Collectors.toList());

        List<HlaSequenceLoci> candidateNucleotideSequences = mRefData.NucleotideSequences.stream()
            .filter(x -> contains(candidateAlleleSpecificProteins, x.getAllele().asFourDigit())).collect(Collectors.toList());

        List<FragmentAlleles> referenceFragmentAlleles = FragmentAlleles.create(
                referenceCoverageFragments, referenceAminoAcidHeterozygousLoci, candidateAminoAcidSequences,
                referenceNucleotideHeterozygousLoci, candidateNucleotideSequences);

        // Complexes
        List<HlaComplex> complexes = HlaComplex.complexes(mConfig, mRefData, referenceFragmentAlleles, candidateAlleles); // , recoveredAlleles, not used

        LL_LOGGER.info("Calculating coverage of {} complexes", complexes.size());
        HlaComplexCoverageFactory coverageFactory = new HlaComplexCoverageFactory(mConfig, mRefData);
        List<HlaComplexCoverage> referenceRankedComplexes = coverageFactory.rankedComplexCoverage(
                mExecutorService, referenceFragmentAlleles, complexes, recoveredAlleles);

        if(referenceRankedComplexes.isEmpty())
        {
            LL_LOGGER.fatal("Failed to calculate complex coverage");
            System.exit(1);
        }

        if(expectedSequences.isEmpty())
        {
            HlaComplexCoverage expectedCoverage = proteinCoverage(
                    referenceFragmentAlleles, expectedSequences.stream().map(x -> x.getAllele()).collect(Collectors.toList()));
            LL_LOGGER.info("Expected allele coverage: {}", expectedCoverage);
        }

        HlaComplexCoverage winningReferenceCoverage = referenceRankedComplexes.get(0).expandToSixAlleles();
        List<HlaAllele> winningAlleles = winningReferenceCoverage.getAlleleCoverage().stream().map(x -> x.Allele).collect(Collectors.toList());
        List<HlaSequenceLoci> winningSequences = candidateSequences.stream()
                .filter(x -> winningAlleles.contains(x.getAllele())).collect(Collectors.toList());

        LL_LOGGER.info("{}", HlaComplexCoverage.header());

        for (HlaComplexCoverage rankedComplex : referenceRankedComplexes)
        {
            LL_LOGGER.info(rankedComplex);
        }

        LL_LOGGER.info("{} - REF - {} CANDIDATES, WINNING ALLELES: {}",
                mConfig.Sample, referenceRankedComplexes.size(),
                HlaAllele.toString(winningReferenceCoverage.getAlleleCoverage().stream().map(x -> x.Allele).collect(Collectors.toList())));

        HlaComplexCoverage winningTumorCoverage = null;
        List<HlaCopyNumber> winningTumorCopyNumber = null;
        List<VariantContextDecorator> somaticVariants = Lists.newArrayList();
        List<SomaticCodingCount> somaticCodingCount = SomaticCodingCount.create(winningAlleles);

        if(!mConfig.TumorBam.isEmpty())
        {
            LL_LOGGER.info("Calculating tumor coverage of winning alleles");

            List<FragmentAlleles> tumorFragmentAlleles = FragmentAlleles.create(aminoAcidPipeline.tumorCoverageFragments(),
                    referenceAminoAcidHeterozygousLoci, candidateAminoAcidSequences,
                    referenceNucleotideHeterozygousLoci, candidateNucleotideSequences);

            winningTumorCoverage = proteinCoverage(tumorFragmentAlleles, winningAlleles).expandToSixAlleles();

            LL_LOGGER.info("Calculating tumor copy number of winning alleles");

            winningTumorCopyNumber = HlaCopyNumber.alleleCopyNumber(winningAlleles, mConfig.GeneCopyNumberFile, winningTumorCoverage);

            // SOMATIC VARIANTS
            somaticVariants.addAll(new SomaticVariants(mConfig).readSomaticVariants());

            if(!somaticVariants.isEmpty())
            {
                LL_LOGGER.info("Calculating somatic variant allele coverage");

                String vcfFilename = mConfig.OutputFilePrefix + ".lilac.somatic.vcf.gz";
                LilacVCF lilacVCF = new LilacVCF(vcfFilename, mConfig.SomaticVcf).writeHeader(version.toString());

                SomaticAlleleCoverage somaticCoverageFactory = new SomaticAlleleCoverage(
                        mConfig, referenceAminoAcidHeterozygousLoci, LOCI_POSITION, somaticVariants, winningSequences);

                for (VariantContextDecorator variant : somaticVariants)
                {
                    List<HlaAlleleCoverage> variantCoverage = somaticCoverageFactory.alleleCoverage(variant, tumorBamReader);
                    Set<HlaAllele> variantAlleles = variantCoverage.stream().map(x -> x.Allele).collect(Collectors.toSet());
                    LL_LOGGER.info("    {} -> {}}", variant, variantCoverage);
                    lilacVCF.writeVariant(variant.context(), variantAlleles);

                    // TODO
                    // somaticCodingCount = somaticCodingCount.addVariant(variant, variantAlleles);
                }
            }
        }
        else
        {
            winningTumorCoverage = HlaComplexCoverage.create(Lists.newArrayList());
            winningTumorCopyNumber = HlaCopyNumber.alleleCopyNumber(winningAlleles);
        }

        HlaOut output = HlaOut.create(winningReferenceCoverage, winningTumorCoverage, winningTumorCopyNumber, somaticCodingCount);

        // QC
        LL_LOGGER.info("Calculating QC Statistics");
        SomaticVariantQC somaticVariantQC = SomaticVariantQC.create(somaticVariants.size(), somaticCodingCount);
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
        candidatesToWrite.add(deflatedSequenceTemplate);

        HlaSequenceLociFile.write(
                String.format("%s.candidates.sequences.txt", mConfig.outputPrefix()),
                A_EXON_BOUNDARIES, B_EXON_BOUNDARIES, C_EXON_BOUNDARIES, candidatesToWrite);

        HlaComplexCoverage.writeToFile(referenceRankedComplexes,
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

    @Override
    public void close()
    {
        mExecutorService.shutdown();
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
