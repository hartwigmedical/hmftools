package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.BachelorProgram.matchesBlacklistExclusion;
import static com.hartwig.hmftools.bachelor.BachelorProgram.matchesWhitelistGeneProtein;
import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONSENSE_OR_FRAMESHIFT;
import static com.hartwig.hmftools.common.variant.CodingEffect.SPLICE;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.CodingEffect;
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

public class ExternalDBFilters
{
    private String mOutputDir;
    private String mInputFilterFile;
    private List<String> mRequiredEffects;
    private List<String> mPanelTranscripts;
    private ProgramBlacklist mConfigBlacklist;
    private ProgramWhitelist mConfigWhitelist;

    private boolean mRunMatching;
    private boolean[] mMatchedBlacklistExclusions;
    private boolean[] mMatchedWhitelistExclusions;

    private BufferedWriter mFilterWriter;

    private static final Logger LOGGER = LogManager.getLogger(ExternalDBFilters.class);

    public ExternalDBFilters()
    {
        mOutputDir = "";
        mInputFilterFile = "";
        mRequiredEffects = Lists.newArrayList();
        mPanelTranscripts = Lists.newArrayList();
        mFilterWriter = null;
        mRunMatching = false;
        mMatchedBlacklistExclusions = null;
        mMatchedWhitelistExclusions = null;
    }

    private static int BACHELOR_FILTER_CSV_FIELD_COUNT = 13;

    public static List<VariantFilter> loadExternalFilters(final String filterFile)
    {
        List<VariantFilter> filters = Lists.newArrayList();

        if (filterFile.isEmpty() || !Files.exists(Paths.get(filterFile)))
            return filters;

        int lineIndex = 0;

        try
        {
            BufferedReader file = new BufferedReader(new FileReader(filterFile));

            file.readLine(); // skip header

            String line = null;

            while ((line = file.readLine()) != null)
            {
                if (line.isEmpty())
                    break;

                ++lineIndex;

                // parse CSV data
                String[] items = line.split(",");

                if (items.length < BACHELOR_FILTER_CSV_FIELD_COUNT)
                {
                    LOGGER.error("invalid item count({}), fileIndex({})", items.length, lineIndex);
                    continue;
                }

                // Gene,TranscriptId,Chromosome,Position,Ref,Alt,CodingEffect,AllEffects,HgvsProtein,HgvsCoding,DBSnpId,ClinvarSignificance,ClinvarSigInfo
                VariantFilter filter = new VariantFilter(
                        items[0],
                        items[1],
                        items[2],
                        Long.parseLong(items[3]),
                        items[4],
                        items[5],
                        CodingEffect.valueOf(items[6]),
                        items[8],
                        items[10],
                        items[11],
                        items[12],
                        -1);

                filters.add(filter);
            }

        }
        catch (IOException e)
        {
            LOGGER.error("Failed to read bachelor input CSV file({}) index({}): {}", filterFile, lineIndex, e.toString());
        }

        LOGGER.info("loaded {} clinvar filter records", filters.size());

        return filters;
    }

    public boolean createFilterFile(final String filterFile, final String outputDir, final Program program)
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

        if(mRunMatching)
        {
            if (!mConfigBlacklist.getExclusion().isEmpty())
            {
                mMatchedBlacklistExclusions = new boolean[mConfigBlacklist.getExclusion().size()];
            }

            if (!mConfigWhitelist.getVariantOrDbSNP().isEmpty())
            {
                mMatchedWhitelistExclusions = new boolean[mConfigWhitelist.getVariantOrDbSNP().size()];
            }
        }

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

                    if(!mMatchedBlacklistExclusions[i])
                    {
                        LOGGER.info("blacklist exclusion {}: {}",
                                mMatchedBlacklistExclusions[i] ? "matched" : "not matched", blacklistExclusionAsString(exclusion));
                    }
                }
            }

            if(mMatchedWhitelistExclusions != null)
            {
                for(int i = 0; i < mMatchedWhitelistExclusions.length; ++i)
                {
                    final Object exclusion = mConfigWhitelist.getVariantOrDbSNP().get(i);

                    if(!mMatchedWhitelistExclusions[i])
                    {
                        LOGGER.info("whitelist exclusion {}: {}",
                                mMatchedWhitelistExclusions[i] ? "matched" : "not matched", whitelistExclusionAsString(exclusion));
                    }
                }
            }
        }
        catch (final TribbleException e)
        {
            LOGGER.error("failed to read filter input VCF file: {}", e.getMessage());
            return false;
        }

        closeBufferedWriter(mFilterWriter);

        return true;
    }

    private static String blacklistExclusionAsString(final ProgramBlacklist.Exclusion blacklist)
    {
        if (blacklist.getHGVSP() != null)
        {
            return String.format("gene(%s) HGVS protein(%s)", blacklist.getGene().getName(), blacklist.getHGVSP().toString());
        }

        if (blacklist.getHGVSC() != null)
        {
            return String.format("gene(%s) HGVS coding(%s)", blacklist.getGene().getName(), blacklist.getHGVSC().toString());
        }

        if(blacklist.getMinCodon() != null)
        {
            return String.format("gene(%s) minCodon(%d)", blacklist.getGene().getName(), blacklist.getMinCodon().intValue());
        }

        if(blacklist.getPosition() != null)
        {
            return String.format("gene(%s) position(%s)", blacklist.getGene().getName(), blacklist.getPosition().toString());
        }

        return "";
    }

    private static String whitelistExclusionAsString(final Object variantOrDbSNP)
    {
        if (variantOrDbSNP instanceof ProgramWhitelist.Variant)
        {
            final ProgramWhitelist.Variant geneProtein = (ProgramWhitelist.Variant) variantOrDbSNP;
            return String.format("gene(%s) protein(%s)", geneProtein.getGene().getName(), geneProtein.getHGVSP().toString());
        }
        else
        {
            return String.format("DbSNP ID(%s)", (String) variantOrDbSNP);
        }
    }


    // Clinvar annotations
    private static String CLINVAR_SIGNIFICANCE = "CLNSIG";
    private static String CLINVAR_SIG_INFO = "CLNSIGCONF";
    private static String CLINVAR_DISEASE_NAME = "CLNDN";
    private static String CLINVAR_MC = "MC";
    private static String CLINVAR_RS_DB_SNP_ID = "RS";
    private static String CLINVAR_PATHOGENIC = "Pathogenic";
    private static String CLINVAR_LIKELY_PATHOGENIC = "Likely_pathogenic";
    private static String CLINVAR_BENIGN = "Benign";
    private static String CLINVAR_LIKELY_BENIGN = "Likely_benign";
    private static String CLINVAR_CONFLICTING = "Conflicting";

    public static boolean isPathogenic(final String clinvarSignificance)
    {
        return clinvarSignificance.contains(CLINVAR_PATHOGENIC) || clinvarSignificance.contains(CLINVAR_LIKELY_PATHOGENIC);
    }

    public static boolean isBenign(final String clinvarSignificance)
    {
        return clinvarSignificance.contains(CLINVAR_BENIGN) || clinvarSignificance.contains(CLINVAR_LIKELY_BENIGN);
    }

    public static boolean isConflicting(final String clinvarSignificance)
    {
        return clinvarSignificance.contains(CLINVAR_CONFLICTING);
    }

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

            String clinvarSignificance = variant.getCommonInfo().getAttributeAsString(CLINVAR_SIGNIFICANCE, "");
            clinvarSignificance = stripArrayChars(clinvarSignificance);

            String clinvarSigInfo = variant.getCommonInfo().getAttributeAsString(CLINVAR_SIG_INFO, "");

            if(!clinvarSigInfo.isEmpty())
            {
                clinvarSigInfo = stripArrayChars(clinvarSigInfo);
            }

            boolean isPathogenic = isPathogenic(clinvarSignificance);

            if(!isPathogenic && isConflicting(clinvarSignificance))
            {
                // look in the significance field for a clear likelihood
                isPathogenic = isPathogenic(clinvarSigInfo) && !isBenign(clinvarSigInfo);
            }

            boolean matchesRequiredEffect = false;

            for (String requiredEffect : mRequiredEffects)
            {
                if (snpEff.effects().contains(requiredEffect))
                {
                    matchesRequiredEffect = true;
                    break;
                }
            }

            if(!matchesRequiredEffect && !isPathogenic)
                continue;

            String gene = snpEff.gene();

            CodingEffect codingEffect = CodingEffect.effect(gene, snpEff.consequences());

            if(codingEffect == NONSENSE_OR_FRAMESHIFT || codingEffect == SPLICE)
            {
                checkExistingBlacklistConditions(gene, variant, snpEff);
            }
            else
            {
                checkExistingWhitelistConditions(gene, variant, snpEff);

                if(!isPathogenic)
                {
                    continue;
                }

                // will form part of the whitelist
            }

            writeFilterRecord(variant, snpEff, gene, codingEffect, clinvarSignificance, clinvarSigInfo);
        }
    }

    private static String stripArrayChars(String str)
    {
        str = str.replaceAll(",", ";");
        str = str.replaceAll("\\[", "");
        str = str.replaceAll("]", "");
        return str;
    }

    private void checkExistingBlacklistConditions(final String gene, final VariantContext variant, final SnpEffAnnotation snpAnnotation)
    {
        if(!mRunMatching)
            return;

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
                    LOGGER.debug("clinar variant for gene({}) matched with blacklist exclusion", gene);
                    mMatchedBlacklistExclusions[index] = true;
                    break;
                }
            }
        }
    }

    private void checkExistingWhitelistConditions(final String gene, final VariantContext variant, final SnpEffAnnotation snpAnnotation)
    {
        if(!mRunMatching)
            return;

        String rsDbSnpId = variant.getCommonInfo().getAttributeAsString(CLINVAR_RS_DB_SNP_ID, "");

        for(int index = 0; index < mConfigWhitelist.getVariantOrDbSNP().size(); ++index)
        {
            final Object variantOrDbSNP = mConfigWhitelist.getVariantOrDbSNP().get(index);

            if (variantOrDbSNP instanceof ProgramWhitelist.Variant)
            {
                final ProgramWhitelist.Variant geneProtein = (ProgramWhitelist.Variant) variantOrDbSNP;

                if (matchesWhitelistGeneProtein(geneProtein, variant, snpAnnotation))
                {
                    LOGGER.debug("clinar variant for gene({}) matched with whitelist geneProtein exclusion", gene);
                    mMatchedWhitelistExclusions[index] = true;
                    break;
                }
            }
            else
            {
                String wlDBSnpId = ((String)variantOrDbSNP).replaceFirst("rs", "");
                if(rsDbSnpId.contains(wlDBSnpId))
                {
                    LOGGER.debug("clinar variant for gene({}) matched with whitelist exclusion rsDbSnpId({})", gene, rsDbSnpId);
                    mMatchedWhitelistExclusions[index] = true;
                    break;
                }
            }
        }
    }

    private void writeFilterRecord(final VariantContext variant, final SnpEffAnnotation snpEff,
            final String gene, CodingEffect codingEffect, final String clinvarSignificance, final String clinvarSigInfo)
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
        String rsDbSnpId = variant.getCommonInfo().getAttributeAsString(CLINVAR_RS_DB_SNP_ID, "");

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
            mFilterWriter.write(String.format("%s,%s,%s,%d,%s,%s,%s,%s",
                    gene, snpEff.transcript(), chromosome, position, ref, alt, codingEffect, effects));

            mFilterWriter.write(String.format(",%s,%s,%s",
                    hgvsp, hgvsc, rsDbSnpId));

            mFilterWriter.write(String.format(",%s,%s,%s,%s",
                    clinvarSignificance, clinvarSigInfo, clinvarDisease, clinvarEffects));

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

            mFilterWriter.write("Gene,TranscriptId,Chromsome,Position,Ref,Alt,CodingEffect,AllEffects");
            mFilterWriter.write(",HgvsProtein,HgvsCoding,DBSnpId");
            mFilterWriter.write(",ClinvarSignificance,ClinvarSigInfo,ClinvarDisease,ClinvarEffects");
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
