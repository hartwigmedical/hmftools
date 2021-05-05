package com.hartwig.hmftools.lilac;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.A_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.LilacConstants.B_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.LilacConstants.C_EXON_BOUNDARIES;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_A;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_B;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_C;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_TRANSCRIPTS;
import static com.hartwig.hmftools.lilac.LilacConstants.LOCI_POSITION;

import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage;
import com.hartwig.hmftools.lilac.coverage.HlaComplexCoverage;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaContext;
import com.hartwig.hmftools.lilac.hla.HlaContextFactory;
import com.hartwig.hmftools.lilac.nuc.NucleotideFragment;
import com.hartwig.hmftools.lilac.nuc.NucleotideFragmentFactory;
import com.hartwig.hmftools.lilac.read.SAMRecordReader;

import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.ThreadFactory;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public class LilacApplication implements AutoCloseable, Runnable
{
    private final long mStartTime;
    private final ThreadFactory mNamedThreadFactory;
    private final ExecutorService mExecutorService;
    // private final NucleotideGeneEnrichment nucleotideGeneEnrichment;
    // private final CommandLine cmd;
    private final LilacConfig mConfig;

    public LilacApplication(final CommandLine cmd, final LilacConfig config)
    {
        mConfig = config;

        mStartTime = System.currentTimeMillis();
        mNamedThreadFactory = null; // new ThreadFactory();;
        mExecutorService = null;
        // NucleotideGeneEnrichment nucleotideGeneEnrichment;

    }

    public void run()
    {
        // load gene definitions
        final Map<String, HmfTranscriptRegion> allTranscripts = HmfGenePanelSupplier.allGenesMap37();
        HLA_TRANSCRIPTS.add(allTranscripts.get(HLA_A));
        HLA_TRANSCRIPTS.add(allTranscripts.get(HLA_B));
        HLA_TRANSCRIPTS.add(allTranscripts.get(HLA_C));
        LOCI_POSITION.initialise(HLA_TRANSCRIPTS);

        String outputQCFile;
        String outputFile;
        List<HlaComplexCoverage> referenceRankedComplexes;
        SequenceCount referenceNucleotideCounts;
        SequenceCount referenceAminoAcidCounts;
        List candidateSequences;
        List expectedSequences;

        LL_LOGGER.info("Starting LILAC with parameters:");
        LL_LOGGER.info("    sample = " + mConfig.Sample);
        LL_LOGGER.info("    minBaseQual = " + mConfig.MinBaseQual);
        LL_LOGGER.info("    minEvidence = " + mConfig.MinEvidence);
        LL_LOGGER.info("    minUniqueCoverage = " + mConfig.MinConfirmedUniqueCoverage);
        LL_LOGGER.info("    minFragmentsPerAllele = " + mConfig.MinFragmentsPerAllele);
        LL_LOGGER.info("    minFragmentsToRemoveSingle = " + mConfig.MinFragmentsToRemoveSingle);

        HlaContextFactory
                hlaContextFactory = new HlaContextFactory(A_EXON_BOUNDARIES, B_EXON_BOUNDARIES, C_EXON_BOUNDARIES);
        HlaContext hlaAContext = hlaContextFactory.hlaA();
        HlaContext hlaBContext = hlaContextFactory.hlaB();
        HlaContext hlaCContext = hlaContextFactory.hlaC();

        ReferenceData refData = new ReferenceData(mConfig.ResourceDir);

        if(!refData.load())
        {
            LL_LOGGER.error("reference data loading failed");
            System.exit(1);
        }

        /*
        LL_LOGGER.info("Reading nucleotide files");
        final List<HlaSequenceLoci> nucleotideSequences = Lists.newArrayList();
        nucleotideSequences.addAll(nucleotideLoci(mConfig.ResourceDir + "/A_nuc.txt"));
        nucleotideSequences.addAll(nucleotideLoci(mConfig.ResourceDir + "/B_nuc.txt"));
        nucleotideSequences.addAll(nucleotideLoci(mConfig.ResourceDir + "/C_nuc.txt"));

        LL_LOGGER.info("Reading protein files");
        final List<HlaSequenceLoci> aminoAcidSequences = Lists.newArrayList();
        aminoAcidSequences.addAll(aminoAcidLoci(mConfig.ResourceDir + "/A_prot.txt"));
        aminoAcidSequences.addAll(aminoAcidLoci(mConfig.ResourceDir + "/B_prot.txt"));
        aminoAcidSequences.addAll(aminoAcidLoci(mConfig.ResourceDir + "/C_prot.txt"));

        final List<HlaSequenceLoci> aminoAcidSequencesWithInserts = aminoAcidSequences.stream().filter(x -> x.containsInserts()).collect(Collectors.toList());
        final List<HlaSequenceLoci> aminoAcidSequencesWithDeletes = aminoAcidSequences.stream().filter(x -> x.containsDeletes()).collect(Collectors.toList());

         */

        LL_LOGGER.info("Querying records from reference bam " + mConfig.ReferenceBam);

        NucleotideFragmentFactory nucleotideFragmentFactory = new NucleotideFragmentFactory(
                mConfig.MinBaseQual, refData.AminoAcidSequencesWithInserts, refData.AminoAcidSequencesWithDeletes, LOCI_POSITION);

        SAMRecordReader tumorBamReader =
                new SAMRecordReader(mConfig.TumorBam, mConfig.RefGenome, HLA_TRANSCRIPTS, nucleotideFragmentFactory);

        SAMRecordReader referenceBamReader =
                new SAMRecordReader(mConfig.ReferenceBam, mConfig.RefGenome, HLA_TRANSCRIPTS, nucleotideFragmentFactory);

        final List<NucleotideFragment> referenceNucleotideFragments = enrichGenes(referenceBamReader.readFromBam());
        final List<NucleotideFragment> tumorNucleotideFragments = Lists.newArrayList();
        if(!mConfig.TumorBam.isEmpty())
        {
            LL_LOGGER.info("Querying records from tumor bam " + mConfig.TumorBam);
            tumorNucleotideFragments.addAll(enrichGenes(tumorBamReader.readFromBam()));
        }

        /*
        amino.AminoAcidFragmentPipeline aminoAcidPipeline =
                new AminoAcidFragmentPipeline(mConfig, referenceNucleotideFragments, (List<? extends nuc.NucleotideFragment>) tumorNucleotideFragments);

        List<amino.AminoAcidFragment> aCandidateFragments = aminoAcidPipeline.referencePhasingFragments(hlaAContext);
        List<amino.AminoAcidFragment> bCandidateFragments = aminoAcidPipeline.referencePhasingFragments(hlaBContext);
        List<amino.AminoAcidFragment> cCandidateFragments = aminoAcidPipeline.referencePhasingFragments(hlaCContext);

        candidates.Candidates candidateFactory = new Candidates(mConfig, nucleotideSequences, aminoAcidSequences);
        List<hla.HlaAllele> aUnphasedCandidates = candidateFactory.unphasedCandidates(hlaAContext, aCandidateFragments);
        List<hla.HlaAllele> bUnphasedCandidates = candidateFactory.unphasedCandidates(hlaBContext, bCandidateFragments);
        List<hla.HlaAllele> cUnphasedCandidates = candidateFactory.unphasedCandidates(hlaCContext, cCandidateFragments);


            List allUnphasedCandidates =
                    CollectionsKt.plus((Collection) CollectionsKt.plus((Collection) aUnphasedCandidates, (Iterable) bUnphasedCandidates), (Iterable) cUnphasedCandidates);
            evidence.PhasedEvidenceFactory phasedEvidenceFactory = new PhasedEvidenceFactory(mConfig);
            List<evidence.PhasedEvidence> aPhasedEvidence = phasedEvidenceFactory.evidence(hlaAContext, aCandidateFragments);
            List<evidence.PhasedEvidence> bPhasedEvidence = phasedEvidenceFactory.evidence(hlaBContext, bCandidateFragments);
            List<PhasedEvidence> cPhasedEvidence = phasedEvidenceFactory.evidence(hlaCContext, cCandidateFragments);
            Iterable iterable2 = $receiver$iv4 = (Iterable) aminoAcidSequences;
            Collection destination$iv$iv3 = new ArrayList();
            for(Object element$iv$iv : $receiver$iv$iv8)
            {
                HlaSequenceLoci it7 = (HlaSequenceLoci) element$iv$iv;
                boolean bl = false;
                if(!mConfig.getExpectedAlleles().contains(it7.getAllele().asFourDigit()))
                {
                    continue;
                }
                destination$iv$iv3.add(element$iv$iv);
            }
            expectedSequences = (List) destination$iv$iv3;
            evidence.PhasedEvidenceValidation.INSTANCE.validateExpected("A", aPhasedEvidence, expectedSequences);
            evidence.PhasedEvidenceValidation.INSTANCE.validateExpected("B", bPhasedEvidence, expectedSequences);
            PhasedEvidenceValidation.INSTANCE.validateExpected("C", cPhasedEvidence, expectedSequences);
            List<hla.HlaAllele> aCandidates = candidateFactory.phasedCandidates(hlaAContext, aUnphasedCandidates, aPhasedEvidence);
            List<hla.HlaAllele> bCandidates = candidateFactory.phasedCandidates(hlaBContext, bUnphasedCandidates, bPhasedEvidence);
            List<hla.HlaAllele> cCandidates = candidateFactory.phasedCandidates(hlaCContext, cUnphasedCandidates, cPhasedEvidence);
            List phasedCandidateAlleles =
                    CollectionsKt.plus((Collection) CollectionsKt.plus((Collection) aCandidates, (Iterable) bCandidates), (Iterable) cCandidates);
            Iterable bl = $receiver$iv3 = (Iterable) mConfig.getStopLossRecoveryAlleles();
            Iterable destination$iv$iv4 = new ArrayList();
            for(Object element$iv$iv : $receiver$iv$iv7)
            {
                it4 = (hla.HlaAllele) element$iv$iv;
                boolean bl2 = false;
                if(!(!phasedCandidateAlleles.contains(it4)))
                {
                    continue;
                }
                destination$iv$iv4.add(element$iv$iv);
            }
            List missingStopLossAlleles = (List) destination$iv$iv4;
            if(referenceBamReader.stopLossOnCIndels() > 0 && !($receiver$iv$iv7 = (Collection) missingStopLossAlleles).isEmpty())
            {
                LL_LOGGER.info("Identified " + referenceBamReader.stopLossOnCIndels() + " stop loss fragments");
                LL_LOGGER.info("    recovered stop loss candidates: "
                        + CollectionsKt.joinToString$default((Iterable) missingStopLossAlleles, (CharSequence) ",", null, null, (int) 0, null, null, (int) 62, null));
                list = missingStopLossAlleles;
            }
            else
            {
                list = CollectionsKt.emptyList();
            }
            List stopLossRecovery = list;
            LL_LOGGER.info("Recovering common un-phased candidates:");
            destination$iv$iv4 = mConfig.getCommonAlleles();
            Iterator $i$f$filter = $receiver$iv2;
            Collection destination$iv$iv5 = new ArrayList();
            for(Object element$iv$iv : $receiver$iv$iv6)
            {
                hla.HlaAllele it8 = (hla.HlaAllele) element$iv$iv;
                boolean bl3 = false;
                if(!(!phasedCandidateAlleles.contains(it8) && allUnphasedCandidates.contains(it8)))
                {
                    continue;
                }
                destination$iv$iv5.add(element$iv$iv);
            }
            List recoveredAlleles = CollectionsKt.plus((Collection) ((List) destination$iv$iv5), (Iterable) stopLossRecovery);
            List candidateAlleles = CollectionsKt.plus((Collection) phasedCandidateAlleles, (Iterable) recoveredAlleles);
            Iterable $receiver$iv7 = aminoAcidSequences;
            it4 = $receiver$iv7;
            Collection destination$iv$iv6 = new ArrayList();
            for(Object element$iv$iv : $receiver$iv$iv5)
            {
                HlaSequenceLoci it9 = (HlaSequenceLoci) element$iv$iv;
                boolean bl4 = false;
                if(!candidateAlleles.contains(it9.getAllele()))
                {
                    continue;
                }
                destination$iv$iv6.add(element$iv$iv);
            }
            candidateSequences = (List) destination$iv$iv6;
            $receiver$iv7 = recoveredAlleles;
            if(!$receiver$iv7.isEmpty())
            {
                LL_LOGGER.info("    recovered " + recoveredAlleles.size() + " common candidate alleles: "
                        + CollectionsKt.joinToString$default((Iterable) recoveredAlleles, (CharSequence) ", ", null, null, (int) 0, null, null, (int) 62, null));
            }
            else
            {
                LL_LOGGER.info("    recovered 0 common candidate alleles");
            }
            List<AminoAcidFragment> referenceCoverageFragments = aminoAcidPipeline.referenceCoverageFragments();
            referenceAminoAcidCounts = SequenceCount.Companion.aminoAcids(minEvidence, referenceCoverageFragments);
            List<Integer> referenceAminoAcidHeterozygousLoci = referenceAminoAcidCounts.heterozygousLoci();
            referenceNucleotideCounts = SequenceCount.Companion.nucleotides(minEvidence, referenceCoverageFragments);
            Set referenceNucleotideHeterozygousLoci =
                    CollectionsKt.intersect((Iterable) referenceNucleotideCounts.heterozygousLoci(), (Iterable) ALL_NUCLEOTIDE_EXON_BOUNDARIES);
            Iterable $i$f$filterTo = $receiver$iv = (Iterable) candidateAlleles;
            Iterable destination$iv$iv7 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
            for(Object item$iv$iv : $receiver$iv$iv4)
            {
                void it10;
                hla.HlaAllele hlaAllele2 = (hla.HlaAllele) item$iv$iv;
                object2 = destination$iv$iv7;
                boolean bl5 = false;
                object = it10.asFourDigit();
                object2.add(object);
            }
            List candidateAlleleSpecificProteins = (List) destination$iv$iv7;
            Iterable $receiver$iv8 = aminoAcidSequences;
            destination$iv$iv7 = $receiver$iv8;
            Iterable destination$iv$iv8 = new ArrayList();
            for(Object element$iv$iv : $receiver$iv$iv3)
            {
                HlaSequenceLoci it11 = (HlaSequenceLoci) element$iv$iv;
                boolean bl6 = false;
                if(!candidateAlleles.contains(it11.getAllele()))
                {
                    continue;
                }
                destination$iv$iv8.add(element$iv$iv);
            }
            List candidateAminoAcidSequences = (List) destination$iv$iv8;
            Iterable $receiver$iv9 = nucleotideSequences;
            destination$iv$iv8 = $receiver$iv9;
            Collection destination$iv$iv9 = new ArrayList();
            for(Object element$iv$iv : $receiver$iv$iv2)
            {
                HlaSequenceLoci it12 = (HlaSequenceLoci) element$iv$iv;
                boolean bl7 = false;
                if(!candidateAlleleSpecificProteins.contains(it12.getAllele().asFourDigit()))
                {
                    continue;
                }
                destination$iv$iv9.add(element$iv$iv);
            }
            List candidateNucleotideSequences = (List) destination$iv$iv9;
            List<read.FragmentAlleles> referenceFragmentAlleles =
                    read.FragmentAlleles.Companion.create(referenceCoverageFragments, (Collection<Integer>) referenceAminoAcidHeterozygousLoci, (Collection<HlaSequenceLoci>) candidateAminoAcidSequences, (Collection<Integer>) referenceNucleotideHeterozygousLoci, (Collection<HlaSequenceLoci>) candidateNucleotideSequences);
            List<coverage.HlaComplex> complexes =
                    HlaComplex.Companion.complexes(mConfig, referenceFragmentAlleles, candidateAlleles, recoveredAlleles);
            LL_LOGGER.info("Calculating coverage of " + complexes.size() + " complexes");
            coverage.HlaComplexCoverageFactory
                    coverageFactory = new coverage.HlaComplexCoverageFactory(mConfig);
            ExecutorService executorService = executorService;
            Intrinsics.checkExpressionValueIsNotNull((Object) executorService, (String) "executorService");
            referenceRankedComplexes =
                    coverageFactory.rankedComplexCoverage(executorService, referenceFragmentAlleles, complexes, recoveredAlleles);
            if(referenceRankedComplexes.isEmpty())
            {
                LL_LOGGER.fatal("Failed to calculate complex coverage");
                int element$iv$iv = 1;
                System.exit(element$iv$iv);
                throw (Throwable) new RuntimeException("System.exit returned normally, while it was supposed to halt JVM.");
            }
            Collection element$iv$iv = expectedSequences;
            if(!element$iv$iv.isEmpty())
            {
                void $receiver$iv$iv11;
                void $receiver$iv10;
                Iterable it12 = expectedSequences;
                object = referenceFragmentAlleles;
                object2 = coverage.HlaComplexCoverageFactory.Companion;
                void bl7 = $receiver$iv10;
                Collection collection2 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv10, (int) 10));
                for(Object item$iv$iv : $receiver$iv$iv11)
                {
                    HlaSequenceLoci hlaSequenceLoci = (HlaSequenceLoci) item$iv$iv;
                    collection = collection2;
                    boolean bl8 = false;
                    hlaAllele = ((HlaSequenceLoci) ((Object) it3)).getAllele();
                    collection.add(hlaAllele);
                }
                collection = (List) collection2;
                coverage.HlaComplexCoverage expectedCoverage =
                        ((coverage.HlaComplexCoverageFactory.Companion) object2).proteinCoverage((List<read.FragmentAlleles>) object, collection);
                LL_LOGGER.info("Expected allele coverage: " + expectedCoverage);
            }
            coverage.HlaComplexCoverage winningReferenceCoverage = referenceRankedComplexes.get(0).expandToSixAlleles();
            List<hla.HlaAllele> winningAlleles = alleles(winningReferenceCoverage.getAlleleCoverage());
            Iterable iterable3 = candidateSequences;
            Iterable $i$f$filter2 = iterable3;
            Collection destination$iv$iv11 = new ArrayList();
            it3 = $receiver$iv$iv.iterator();
            while(it3.hasNext())
            {
                Object element$iv$iv2 = it3.next();
                candidate = (HlaSequenceLoci) element$iv$iv2;
                boolean bl9 = false;
                if(!winningAlleles.contains(((HlaSequenceLoci) candidate).getAllele()))
                {
                    continue;
                }
                destination$iv$iv11.add(element$iv$iv2);
            }
            Set winningSequences = CollectionsKt.toSet((Iterable) ((List) destination$iv$iv11));
            LL_LOGGER.info(coverage.HlaComplexCoverage.Companion.header());
            for(coverage.HlaComplexCoverage hlaComplexCoverage : referenceRankedComplexes)
            {
                LL_LOGGER.info((Object) hlaComplexCoverage);
            }
            Iterable iterable4 = winningReferenceCoverage.getAlleleCoverage();
            object = new StringBuilder().append(mConfig.getSample())
                    .append(" - REF - ")
                    .append(referenceRankedComplexes.size())
                    .append(" CANDIDATES, WINNING ALLELES: ");
            object2 = LL_LOGGER;
            $receiver$iv$iv = iterable4;
            Collection destination$iv$iv10 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) iterable4, (int) 10));
            it3 = $receiver$iv$iv.iterator();
            while(it3.hasNext())
            {
                Object item$iv$iv = it3.next();
                candidate = (coverage.HlaAlleleCoverage) item$iv$iv;
                collection = destination$iv$iv10;
                boolean bl10 = false;
                hlaAllele = ((coverage.HlaAlleleCoverage) ((Object) it2)).getAllele();
                collection.add(hlaAllele);
            }
            collection = (List) destination$iv$iv10;
            object2.info(((StringBuilder) object).append(collection).toString());
            Object var58_83 = null;
            List<cna.HlaCopyNumber> winningTumorCopyNumber = null;
            List<VariantContextDecorator> somaticVariants = null;
            List<variant.SomaticCodingCount> somaticCodingCount = variant.SomaticCodingCount.Companion.create(winningAlleles);
            CharSequence item$iv$iv = mConfig.getTumorBam();
            if(item$iv$iv.length() > 0)
            {
                LL_LOGGER.info("Calculating tumor coverage of winning alleles");
                List<read.FragmentAlleles> tumorFragmentAlleles =
                        FragmentAlleles.Companion.create(aminoAcidPipeline.tumorCoverageFragments(), (Collection<Integer>) referenceAminoAcidHeterozygousLoci, (Collection<HlaSequenceLoci>) candidateAminoAcidSequences, (Collection<Integer>) referenceNucleotideHeterozygousLoci, (Collection<HlaSequenceLoci>) candidateNucleotideSequences);
                coverage.HlaComplexCoverage hlaComplexCoverage =
                        HlaComplexCoverageFactory.Companion.proteinCoverage(tumorFragmentAlleles, (Collection<hla.HlaAllele>) winningAlleles)
                                .expandToSixAlleles();
                LL_LOGGER.info("Calculating tumor copy number of winning alleles");
                winningTumorCopyNumber =
                        cna.HlaCopyNumber.Companion.alleleCopyNumber(winningAlleles, mConfig.getGeneCopyNumberFile(), hlaComplexCoverage);

                somaticVariants = new SomaticVariants(mConfig).readSomaticVariants();

                it2 = somaticVariants;
                if(!it2.isEmpty())
                {
                    LL_LOGGER.info("Calculating somatic variant allele coverage");

                    variant.LilacVCF lilacVCF = new LilacVCF(
                            mConfig.OutputFilePrefix + ".lilac.somatic.vcf.gz", mConfig.SomaticVcf).writeHeader();

                    variant.SomaticAlleleCoverage somaticCoverageFactory =
                            new SomaticAlleleCoverage(mConfig, (Collection<Integer>) referenceAminoAcidHeterozygousLoci, LOCI_POSITION, somaticVariants, winningSequences);

                    for(VariantContextDecorator variant : somaticVariants)
                    {
                        List<coverage.HlaAlleleCoverage> variantCoverage = somaticCoverageFactory.alleleCoverage(variant, tumorBamReader);
                        Set<hla.HlaAllele> variantAlleles = Sets.newHashSet();
                        variantCoverage.forEach(x -> variantAlleles.add(x.getAllele());

                        LL_LOGGER.info("    " + variant + " -> " + variantCoverage);
                        VariantContext variantContext = variant.context();
                        lilacVCF.writeVariant(variantContext, variantAlleles);

                        somaticCodingCount = SomaticCodingCount.Companion.addVariant(somaticCodingCount, variant, variantAlleles);
                    }
                }
            }
            else
            {
                coverage.HlaComplexCoverage
                        hlaComplexCoverage = coverage.HlaComplexCoverage.Companion.create(CollectionsKt.emptyList());
                winningTumorCopyNumber = HlaCopyNumber.Companion.alleleCopyNumber(winningAlleles);
                somaticVariants = CollectionsKt.emptyList();
            }
            out.HlaOut output =
                    HlaOut.Companion.create(winningReferenceCoverage, (coverage.HlaComplexCoverage) var58_86, winningTumorCopyNumber, somaticCodingCount);
            LL_LOGGER.info("Calculating QC Statistics");
            qc.SomaticVariantQC
                    somaticVariantQC = SomaticVariantQC.Companion.create(somaticVariants.size(), somaticCodingCount);
            qc.AminoAcidQC aminoAcidQC = AminoAcidQC.Companion.create(winningSequences, referenceAminoAcidCounts);
            qc.HaplotypeQC haplotypeQC =
                    HaplotypeQC.Companion.create(3, winningSequences, CollectionsKt.plus((Collection) CollectionsKt.plus((Collection) aPhasedEvidence, (Iterable) bPhasedEvidence), (Iterable) cPhasedEvidence), referenceAminoAcidCounts);
            qc.BamQC bamQC = BamQC.Companion.create(referenceBamReader);
            CoverageQC coverageQC = CoverageQC.Companion.create(referenceNucleotideFragments.size(), winningReferenceCoverage);
            qc.LilacQC
                    lilacQC = LilacQC.Companion.create(aminoAcidQC, bamQC, coverageQC, haplotypeQC, somaticVariantQC);
            LL_LOGGER.info("QC Stats:");
            LL_LOGGER.info("    "
                    + CollectionsKt.joinToString$default((Iterable) lilacQC.header(), (CharSequence) ",", null, null, (int) 0, null, null, (int) 62, null));
            LL_LOGGER.info("    "
                    + CollectionsKt.joinToString$default((Iterable) lilacQC.body(), (CharSequence) ",", null, null, (int) 0, null, null, (int) 62, null));
            LL_LOGGER.info("Writing output to " + outputDir);
            outputFile = mConfig.getOutputFilePrefix() + ".lilac.txt";
            outputQCFile = mConfig.getOutputFilePrefix() + ".lilac.qc.txt";
            output.write(outputFile);
            lilacQC.writefile(outputQCFile);
            Iterable $receiver$iv12 = aminoAcidSequences;
            for(Object element$iv22 : $receiver$iv12)
            {
                it = (HlaSequenceLoci) element$iv22;
                boolean bl11 = false;
                if(!Intrinsics.areEqual((Object) ((HlaSequenceLoci) it).getAllele(), (Object) DEFLATE_TEMPLATE))
                {
                    continue;
                }
                break block27;
            }
            throw (Throwable) new NoSuchElementException("Collection contains no element matching the predicate.");
        }
        HlaSequenceLoci deflatedSequenceTemplate = (HlaSequenceLoci) element$iv22;
        Iterable $receiver$iv =
                CollectionsKt.distinct((Iterable) CollectionsKt.plus((Collection) CollectionsKt.plus((Collection) candidateSequences, (Iterable) expectedSequences), (Object) deflatedSequenceTemplate));
        element$iv22 = $receiver$iv;
        it = new Comparator<T>()
        {

            public final int compare(T a, T b)
            {
                HlaSequenceLoci it = (HlaSequenceLoci) a;
                boolean bl = false;
                Comparable comparable = it.getAllele();
                it = (seq.HlaSequenceLoci) b;
                Comparable comparable2 = comparable;
                bl = false;
                hla.HlaAllele hlaAllele = it.getAllele();
                return ComparisonsKt.compareValues((Comparable) comparable2, (Comparable) hlaAllele);
            }
        };
        List candidateToWrite = CollectionsKt.sortedWith(element$iv22, (Comparator) it);
        HlaSequenceLociFile.INSTANCE.write(outputDir + '/' + sample
                + ".candidates.sequences.txt", A_EXON_BOUNDARIES, B_EXON_BOUNDARIES, C_EXON_BOUNDARIES, candidateToWrite);
        HlaComplexCoverage.Companion.writeToFile(referenceRankedComplexes, outputDir + '/' + sample + ".candidates.coverage.txt");
        referenceAminoAcidCounts.writeVertically(outputDir + '/' + sample + ".candidates.aminoacids.txt");
        referenceNucleotideCounts.writeVertically(outputDir + '/' + sample + ".candidates.nucleotides.txt");
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

    private final List<HlaAllele> alleles(final List<HlaAlleleCoverage> minConfirmedUniqueCoverage)
    {
        return Lists.newArrayList();

        /*
        void $receiver$iv$iv;
        Iterable $receiver$iv;
        Iterable iterable = $receiver$iv = (Iterable) minConfirmedUniqueCoverage;
        Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            void it;
            coverage.HlaAlleleCoverage hlaAlleleCoverage = (HlaAlleleCoverage) item$iv$iv;
            Collection collection = destination$iv$iv;
            boolean bl = false;
            HlaAllele hlaAllele = it.getAllele();
            collection.add(hlaAllele);
        }
        return (List) destination$iv$iv;

         */
    }

    private final List<NucleotideFragment> enrichGenes(final List<NucleotideFragment> minConfirmedUniqueCoverage)
    {
        return Lists.newArrayList();
        // TODO
        // return nucleotideGeneEnrichment.enrich(minConfirmedUniqueCoverage);
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final VersionInfo version = new VersionInfo("lilac.version");
        LL_LOGGER.info("Lilac version: {}", version.version());

        final Options options = LilacConfig.createOptions();
        final CommandLine cmd = createCommandLine(args, options);

        // if(cmd.hasOption(LOG_DEBUG))
        //    Configurator.setRootLevel(Level.DEBUG);

        LilacApplication lilac = new LilacApplication(cmd, new LilacConfig(cmd));
        lilac.run();

        /*
            } catch (e: IOException) {
            logger.warn(e)
            exitProcess(1)
        } catch (e: ParseException) {
            logger.warn(e)
            val formatter = HelpFormatter()
            formatter.printHelp("lilac", options)
            exitProcess(1)
        }

         */
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
