package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.EligibilityReport.MatchType.GENE_TRANSCRIPT;
import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.VariantConsequence.FRAMESHIFT_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SPLICE_ACCEPTOR_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.SPLICE_DONOR_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.STOP_GAINED;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory;

import nl.hartwigmedicalfoundation.bachelor.GeneIdentifier;
import nl.hartwigmedicalfoundation.bachelor.Program;
import nl.hartwigmedicalfoundation.bachelor.ProgramPanel;
import nl.hartwigmedicalfoundation.bachelor.SnpEffect;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.TribbleException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class FilterFileBuilder
{
    private String mOutputDir;
    private String mInputFilterFile;
    private List<String> mRequiredEffects;
    private List<String> mPanelTranscripts;
    private BufferedWriter mFilterWriter;

    private static final Logger LOGGER = LogManager.getLogger(FilterFileBuilder.class);

    public FilterFileBuilder()
    {
        mOutputDir = "";
        mInputFilterFile = "";
        mRequiredEffects = Lists.newArrayList();
        mPanelTranscripts = Lists.newArrayList();
        mFilterWriter = null;
    }

    public boolean initialise(final String filterFile, final String outputDir, final Program program)
    {
        mOutputDir = outputDir;
        mInputFilterFile = filterFile;

        if (!Files.exists(Paths.get(mInputFilterFile)))
        {
            LOGGER.error("failed to load filter input file: {}", mInputFilterFile);
            return false;
        }

        if(!initialiseFilterWriter())
            return false;

        ProgramPanel programConfig = program.getPanel().get(0);

        List<GeneIdentifier> geneslist = programConfig.getGene();

        // take up a collection of the effects to search for
        mRequiredEffects = programConfig.getSnpEffect().stream().map(SnpEffect::value).collect(Collectors.toList());
        mPanelTranscripts = geneslist.stream().map(GeneIdentifier::getEnsembl).collect(Collectors.toList());

        return true;
    }

    public void run()
    {
        try
        {
            File vcfFile = new File(mInputFilterFile);

            final VCFFileReader vcfReader = new VCFFileReader(vcfFile, false);

            for (final VariantContext variant : vcfReader)
            {
                processVariant(variant);
            }
        }
        catch (final TribbleException e)
        {
            LOGGER.error("failed to read filter input VCF file: {}", e.getMessage());
            return;
        }

        closeBufferedWriter(mFilterWriter);
    }

    // Clinvar annotations
    private static String CLINVAR_SIGNIFICANCE = "CLNSIG";
    private static String CLINVAR_DISEASE_NAME = "CLNDN";
    private static String CLINVAR_MC = "MC";
    private static String CLINVAR_PATHOGENIC = "Pathogenic";
    private static String CLINVAR_LIKELY_PATHOGENIC = "Likely_pathogenic";
    private static String CLINVAR_BENIGN = "Benign";
    private static String CLINVAR_LIKELY_BENIGN = "Likely_benign";

    private void processVariant(final VariantContext variant)
    {
        // LOGGER.debug("read var({}) chr({}))", variant.getID(), variant.getContig());

        // first check against the genes list
        final List<SnpEffAnnotation> sampleAnnotations = SnpEffAnnotationFactory.fromContext(variant);

        for (int i = 0; i < sampleAnnotations.size(); ++i)
        {
            final SnpEffAnnotation snpEff = sampleAnnotations.get(i);

            if (!snpEff.isTranscriptFeature())
                continue;

            if (!mPanelTranscripts.contains(snpEff.transcript()))
                continue;

            boolean matchesRequiredEffect = false;

            for (String requiredEffect : mRequiredEffects)
            {
                if (snpEff.effects().contains(requiredEffect))
                {
                    matchesRequiredEffect = true;
                    break;
                }
            }

            if(!matchesRequiredEffect)
                continue;

            // next check whether the significance makes it irrelevant
            String effects = snpEff.effects();
            List<String> effectsList = Arrays.stream(effects.split("&")).collect(Collectors.toList());
            String clinvarSignificance = variant.getCommonInfo().getAttributeAsString(CLINVAR_SIGNIFICANCE, "");

            boolean hasCodingEffect = hasTypicalCodingEffect(effectsList);

            if(hasCodingEffect)
            {
                if(clinvarSignificance.equals(CLINVAR_BENIGN) || clinvarSignificance.equals(CLINVAR_LIKELY_BENIGN))
                {
                    continue;
                }
            }
            else
            {
                if(!clinvarSignificance.equals(CLINVAR_PATHOGENIC) && !clinvarSignificance.equals(CLINVAR_LIKELY_PATHOGENIC))
                {
                    continue;
                }
            }

            String gene = snpEff.gene();
            String transcriptId = snpEff.transcript();

            String chromosome = variant.getContig();
            long position = variant.getStart();
            String ref = variant.getReference().getBaseString();
            String alt = variant.getAlleles().get(1).getBaseString();

            ref = ref.replaceAll("\\*", "");
            alt = alt.replaceAll("\\*", "");

            String clinvarDisease = variant.getCommonInfo().getAttributeAsString(CLINVAR_DISEASE_NAME, "");
            String clinvarEffects = variant.getCommonInfo().getAttributeAsString(CLINVAR_MC, "");

            if(LOGGER.isDebugEnabled())
            {
                // now extract other required Clinvar info
                LOGGER.debug("var({}:{}) ref({}) alt({}) effect({}) gene({} trans={}) clinvar({}, {}, {})",
                        variant.getContig(), variant.getStart(),
                        variant.getReference().getBaseString(), variant.getAlleles().get(1).getBaseString(),
                        snpEff.effects(), snpEff.transcript(),
                        clinvarSignificance,  clinvarDisease, clinvarEffects);
            }

            writeFilterRecord(
                    chromosome, position, ref, alt, gene, transcriptId,
                    effects, clinvarDisease, clinvarSignificance, clinvarEffects);
        }
    }

    private boolean hasTypicalCodingEffect(List<String> effectsList)
    {
        for(final String effect : effectsList)
        {
            if(FRAMESHIFT_VARIANT.isParentTypeOf(effect) || STOP_GAINED.isParentTypeOf(effect)
            || SPLICE_ACCEPTOR_VARIANT.isParentTypeOf(effect) || SPLICE_DONOR_VARIANT.isParentTypeOf(effect))
            {
                return true;
            }
        }

        return false;
    }

    private void writeFilterRecord(final String chromosome, long position, final String ref, final String alt, final String gene,
            final String transcriptId, final String effect, final String clinvarDisease, final String clinvarSignificance,
            final String clinvarEffects)
    {
        try
        {
            mFilterWriter.write(String.format("%s,%s,%s,%d,%s,%s,%s",
                    gene, transcriptId, chromosome, position, ref, alt, effect, ""));

            mFilterWriter.write(String.format(",%s,%s,%s",
                    clinvarSignificance, clinvarDisease, clinvarEffects));
            mFilterWriter.newLine();

        }
        catch(IOException e)
        {

        }
    }

    private boolean initialiseFilterWriter()
    {
        try
        {
            String outputDir = mOutputDir;

            if (!outputDir.endsWith(File.separator))
                outputDir += File.separator;

            String filterFileName = outputDir + "BACHELOR_CLINVAR_FILTERS.csv";

            mFilterWriter = createBufferedWriter(filterFileName, false);

            mFilterWriter.write("Gene,TranscriptId,Chromsome,Position,Ref,Alt,Effect,ProteinCodon");
            mFilterWriter.write(",ClinvarDiagnosis,ClinvarDisesase,ClinvarEffects");
            mFilterWriter.newLine();
        }
        catch (IOException e)
        {
            LOGGER.error("failed to create output filter file: {}", e.toString());
            return false;
        }

        return true;
    }


}
