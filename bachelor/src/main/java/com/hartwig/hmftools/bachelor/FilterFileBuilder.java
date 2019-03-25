package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.EligibilityReport.MatchType.GENE_TRANSCRIPT;
import static com.hartwig.hmftools.bachelor.predicates.BlacklistPredicate.asString;
import static com.hartwig.hmftools.bachelor.predicates.BlacklistPredicate.matchesBlacklistExclusion;
import static com.hartwig.hmftools.bachelor.predicates.BlacklistPredicate.proteinPosition;
import static com.hartwig.hmftools.bachelor.predicates.WhitelistPredicate.matchesWhitelistDbSNPId;
import static com.hartwig.hmftools.bachelor.predicates.WhitelistPredicate.matchesWhitelistGeneProtein;
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

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bachelor.predicates.WhitelistPredicate;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory;

import nl.hartwigmedicalfoundation.bachelor.GeneIdentifier;
import nl.hartwigmedicalfoundation.bachelor.Program;
import nl.hartwigmedicalfoundation.bachelor.ProgramBlacklist;
import nl.hartwigmedicalfoundation.bachelor.ProgramPanel;
import nl.hartwigmedicalfoundation.bachelor.ProgramWhitelist;
import nl.hartwigmedicalfoundation.bachelor.SnpEffect;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.tribble.TribbleException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class FilterFileBuilder
{
    private String mOutputDir;
    private String mInputFilterFile;
    private List<String> mRequiredEffects;
    private List<String> mPanelTranscripts;
    private ProgramBlacklist mConfigBlacklist;
    private ProgramWhitelist mConfigWhitelist;
    private boolean[] mMatchedBlacklistExclusions;
    private boolean[] mMatchedWhitelistExclusions;

    private BufferedWriter mFilterWriter;

    private static final Logger LOGGER = LogManager.getLogger(FilterFileBuilder.class);

    public FilterFileBuilder()
    {
        mOutputDir = "";
        mInputFilterFile = "";
        mRequiredEffects = Lists.newArrayList();
        mPanelTranscripts = Lists.newArrayList();
        mFilterWriter = null;
        mMatchedBlacklistExclusions = null;
        mMatchedWhitelistExclusions = null;
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
        mConfigBlacklist = program.getBlacklist();
        mConfigWhitelist = program.getWhitelist();

        if(!mConfigBlacklist.getExclusion().isEmpty())
        {
            mMatchedBlacklistExclusions = new boolean[mConfigBlacklist.getExclusion().size()];
        }

        if(!mConfigWhitelist.getVariantOrDbSNP().isEmpty())
        {
            mMatchedWhitelistExclusions = new boolean[mConfigWhitelist.getVariantOrDbSNP().size()];
        }

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

            // log any unmap'tched blacklist or whitelist items
            if(mMatchedBlacklistExclusions != null)
            {
                for(int i = 0; i < mMatchedBlacklistExclusions.length; ++i)
                {
                    final ProgramBlacklist.Exclusion exclusion = mConfigBlacklist.getExclusion().get(i);

                    LOGGER.info("blacklist exclusion {}: {}",
                            mMatchedBlacklistExclusions[i] ? "matched" : "not matched", asString(exclusion));
                }
            }

            if(mMatchedWhitelistExclusions != null)
            {
                for(int i = 0; i < mMatchedWhitelistExclusions.length; ++i)
                {
                    final Object exclusion = mConfigWhitelist.getVariantOrDbSNP().get(i);

                    LOGGER.info("whitelist exclusion {}: {}",
                            mMatchedWhitelistExclusions[i] ? "matched" : "not matched",
                            WhitelistPredicate.asString(exclusion));
                }
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

            String gene = snpEff.gene();

            // next check whether the significance makes it irrelevant
            String effects = snpEff.effects();
            List<String> effectsList = Arrays.stream(effects.split("&")).collect(Collectors.toList());
            String clinvarSignificance = variant.getCommonInfo().getAttributeAsString(CLINVAR_SIGNIFICANCE, "");

            boolean hasCodingEffect = hasTypicalCodingEffect(effectsList);

            if(hasCodingEffect)
            {
                checkExistingBlacklistConditions(gene, variant, snpEff);

                if(clinvarSignificance.equals(CLINVAR_BENIGN) || clinvarSignificance.equals(CLINVAR_LIKELY_BENIGN))
                {
                    continue;
                }
            }
            else
            {
                checkExistingWhitelistConditions(gene, variant, snpEff);

                if(!clinvarSignificance.equals(CLINVAR_PATHOGENIC) && !clinvarSignificance.equals(CLINVAR_LIKELY_PATHOGENIC))
                {
                    continue;
                }
            }

            writeFilterRecord(variant, snpEff, gene, clinvarSignificance);
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

    private boolean checkExistingBlacklistConditions(final String gene, final VariantContext variant, final SnpEffAnnotation snpAnnotation)
    {
        // check whether any of the configured blacklist conditions are covered by the clinvar variants
        /*
         <Exclusion><Gene name="BRCA2"/>
           <MinCodon>3326</MinCodon> OR <HGVS.c>1094_1095insAATT</HGVS.c> OR <Position>11:108121410</Position>
         </Exclusion>
        */

        for(int index = 0; index < mConfigBlacklist.getExclusion().size(); ++index)
        {
            final ProgramBlacklist.Exclusion exclusion = mConfigBlacklist.getExclusion().get(index);

            if(exclusion.getGene().getName().equals(gene))
            {
                if(matchesBlacklistExclusion(exclusion, variant, snpAnnotation))
                {
                    /*
                    List<Integer> proteinPositions = proteinPosition(snpAnnotation);
                    String hgvsp = snpAnnotation.hgvsProtein();
                    String hgvsc = snpAnnotation.hgvsCoding();
                    */

                    LOGGER.debug("clinar variant for gene({}) matched with blacklist exclusion", gene);
                    mMatchedBlacklistExclusions[index] = true;
                    return true;
                }
            }
        }

        return false;
    }

    private boolean checkExistingWhitelistConditions(final String gene, final VariantContext variant, final SnpEffAnnotation snpAnnotation)
    {
        for(int index = 0; index < mConfigWhitelist.getVariantOrDbSNP().size(); ++index)
        {
            final Object variantOrDbSNP = mConfigWhitelist.getVariantOrDbSNP().get(index);

            if (variantOrDbSNP instanceof ProgramWhitelist.Variant)
            {
                final ProgramWhitelist.Variant geneProtein = (ProgramWhitelist.Variant) variantOrDbSNP;

                if (matchesWhitelistGeneProtein(geneProtein, variant, snpAnnotation))
                {
                    LOGGER.debug("clinar variant for gene({}) matched with whitelist exclusion", gene);
                    mMatchedWhitelistExclusions[index] = true;
                    return true;
                }
            }
            else
            {
                if(matchesWhitelistDbSNPId((String)variantOrDbSNP, variant, snpAnnotation))
                {
                    LOGGER.debug("clinar variant for gene({}) matched with whitelist exclusion", gene);
                    mMatchedWhitelistExclusions[index] = true;
                    return true;
                }
            }
        }

        return false;
    }

    private void writeFilterRecord(final VariantContext variant, final SnpEffAnnotation snpEff,
            final String gene, final String clinvarSignificance)
    {
        String transcriptId = snpEff.transcript();
        String chromosome = variant.getContig();
        long position = variant.getStart();
        String ref = variant.getReference().getBaseString();
        String alt = variant.getAlleles().get(1).getBaseString();

        ref = ref.replaceAll("\\*", "");
        alt = alt.replaceAll("\\*", "");

        String effects = snpEff.effects();
        String clinvarDisease = variant.getCommonInfo().getAttributeAsString(CLINVAR_DISEASE_NAME, "");
        String clinvarEffects = variant.getCommonInfo().getAttributeAsString(CLINVAR_MC, "");

        String hgvsp = snpEff.hgvsProtein();
        String hgvsc = snpEff.hgvsCoding();

        if(LOGGER.isDebugEnabled())
        {
            // now extract other required Clinvar info
            LOGGER.debug("var({}:{}) ref({}) alt({}) effect({}) gene({} trans={}) clinvar({}, {}, {})",
                    variant.getContig(), position, ref, alt, snpEff.effects(), gene, transcriptId,
                    clinvarSignificance,  clinvarDisease, clinvarEffects);
        }

        try
        {
            mFilterWriter.write(String.format("%s,%s,%s,%d,%s,%s,%s",
                    gene, snpEff.transcript(), chromosome, position, ref, alt, effects, ""));

            mFilterWriter.write(String.format(",%s,%s,%s,%s,%s",
                    hgvsp, hgvsc, clinvarSignificance, clinvarDisease, clinvarEffects));
            mFilterWriter.newLine();

        }
        catch(IOException e)
        {
            LOGGER.error("error writing filter output: {}", e.toString());
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
            mFilterWriter.write(",HgvsProtein,HgvsCoding,ClinvarDiagnosis,ClinvarDisesase,ClinvarEffects");
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
