package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.BachelorApplication.BATCH_OUTPUT_DIR;
import static com.hartwig.hmftools.bachelor.BachelorApplication.CONFIG_XML;
import static com.hartwig.hmftools.bachelor.BachelorApplication.LOG_DEBUG;
import static com.hartwig.hmftools.bachelor.BachelorApplication.createCommandLine;
import static com.hartwig.hmftools.bachelor.BachelorApplication.loadXML;
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
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bachelor.datamodel.GeneIdentifier;
import com.hartwig.hmftools.bachelor.datamodel.Program;
import com.hartwig.hmftools.bachelor.datamodel.ProgramPanel;
import com.hartwig.hmftools.bachelor.datamodel.SnpEffect;
import com.hartwig.hmftools.bachelor.types.VariantFilter;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;

import htsjdk.tribble.TribbleException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class ExternalDBFilters
{
    private String mInputFilterFile;
    private List<String> mRequiredEffects;
    private List<String> mPanelTranscripts;

    private BufferedWriter mFilterWriter;

    private static final String OUTPUT_DIR = "output_dir";
    private static final String CREATE_FILTER_FILE = "create_filter_file";

    private static final Logger LOGGER = LogManager.getLogger(ExternalDBFilters.class);

    public static void main(final String... args)
    {
        final Options options = new Options();
        options.addOption(CONFIG_XML, true, "XML with genes, black and white lists");
        options.addOption(OUTPUT_DIR, true, "Optional: when in batch mode, all output written to single file");
        options.addOption(CREATE_FILTER_FILE, true, "Optional: create black and white list filter files");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");

        try
        {
            final CommandLine cmd = createCommandLine(options, args);

            if (cmd.hasOption(LOG_DEBUG))
                Configurator.setRootLevel(Level.DEBUG);

            LOGGER.info("building Clinvar filter files");
            final String filterInputFile = cmd.getOptionValue(CREATE_FILTER_FILE);

            String outputDir = cmd.getOptionValue(OUTPUT_DIR, "");

            if (!cmd.hasOption(CONFIG_XML))
                return;

            Map<String, Program> configMap = null;

            try
            {
                configMap = loadXML(Paths.get(cmd.getOptionValue(CONFIG_XML)));
            }
            catch(Exception e)
            {
                LOGGER.error("error loading XML: {}", e.toString());
                return;
            }

            final Program program = configMap.values().iterator().next();

            ExternalDBFilters filterFileBuilder = new ExternalDBFilters(filterInputFile);
            filterFileBuilder.createFilterFile(program, outputDir);

            LOGGER.info("filter file creation complete");
        }
        catch(ParseException e)
        {
            LOGGER.error("config error: {}", e.toString());
            return;
        }
    }


    public ExternalDBFilters(final String filterInputFile)
    {
        mInputFilterFile = filterInputFile;
        mRequiredEffects = Lists.newArrayList();
        mPanelTranscripts = Lists.newArrayList();
        mFilterWriter = null;
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

    public boolean createFilterFile(final Program program, final String outputDir)
    {
        if (!Files.exists(Paths.get(mInputFilterFile)))
        {
            LOGGER.error("failed to load filter input file: {}", mInputFilterFile);
            return false;
        }

        if(!initialiseFilterWriter(outputDir))
            return false;

        ProgramPanel programConfig = program.getPanel().get(0);

        List<GeneIdentifier> geneslist = programConfig.getGene();

        // take up a collection of the effects to search for
        mRequiredEffects = programConfig.getSnpEffect().stream().map(SnpEffect::value).collect(Collectors.toList());
        mPanelTranscripts = geneslist.stream().map(GeneIdentifier::getEnsembl).collect(Collectors.toList());

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
            return false;
        }

        closeBufferedWriter(mFilterWriter);

        return true;
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
                //checkExistingBlacklistConditions(gene, variant, snpEff);
            }
            else
            {
                //checkExistingWhitelistConditions(gene, variant, snpEff);

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

    private boolean initialiseFilterWriter(final String outputDir)
    {
        try
        {
            String filterFileName = outputDir;

            if (!filterFileName.endsWith(File.separator))
                filterFileName += File.separator;

            filterFileName += "BACHELOR_CLINVAR_FILTERS.csv";

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
