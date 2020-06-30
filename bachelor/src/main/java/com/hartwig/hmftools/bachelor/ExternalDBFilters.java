package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.types.BachelorConfig.CONFIG_XML;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.LOG_DEBUG;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.loadXML;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
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
import com.google.common.collect.Maps;
import com.hartwig.hmftools.bachelor.datamodel.GeneIdentifier;
import com.hartwig.hmftools.bachelor.datamodel.Program;
import com.hartwig.hmftools.bachelor.datamodel.ProgramPanel;
import com.hartwig.hmftools.bachelor.datamodel.SnpEffect;
import com.hartwig.hmftools.bachelor.types.VariantFilter;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationFactory;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;

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

    public static void main(final String... args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(CONFIG_XML, true, "XML with genes, black and white lists");
        options.addOption(OUTPUT_DIR, true, "Optional, directory to where the output ClinVar filter file will be written to.");
        options.addOption(CREATE_FILTER_FILE, true, "Input file to create the ClinVar db filters");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");

        final CommandLine cmd = new DefaultParser().parse(options, args);

        if (!cmd.hasOption(CONFIG_XML))
        {
            LOGGER.error("mising XML config file");
            return;
        }

        if (cmd.hasOption(LOG_DEBUG))
            Configurator.setRootLevel(Level.DEBUG);

        LOGGER.info("generating ClinVar filter file");

        Map<String, Program> configMap = Maps.newHashMap();
        loadXML(Paths.get(cmd.getOptionValue(CONFIG_XML)), configMap);

        Program program = configMap.values().iterator().next();

        String filterInputFile = cmd.getOptionValue(CREATE_FILTER_FILE);
        ExternalDBFilters filterFileBuilder = new ExternalDBFilters(filterInputFile);

        String outputDir = cmd.getOptionValue(OUTPUT_DIR, "");
        filterFileBuilder.createFilterFile(program, outputDir);

        LOGGER.info("filter file creation complete");
    }

    private ExternalDBFilters(String filterInputFile)
    {
        mInputFilterFile = filterInputFile;
        mRequiredEffects = Lists.newArrayList();
        mPanelTranscripts = Lists.newArrayList();
        mFilterWriter = null;
    }

    private static final int BACHELOR_FILTER_CSV_FIELD_COUNT = 13;

    static List<VariantFilter> loadExternalFilters(String filterFile)
    {
        List<VariantFilter> filters = Lists.newArrayList();

        if (filterFile.isEmpty() || !Files.exists(Paths.get(filterFile)))
            return filters;

        int lineIndex = 0;

        try
        {
            BufferedReader file = new BufferedReader(new FileReader(filterFile));

            final String headers = file.readLine();
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(headers, ",");

            int geneIndex = fieldsIndexMap.get("Gene");
            int transcriptIndex = fieldsIndexMap.get("TranscriptId");
            int chromosomeIndex = fieldsIndexMap.get("Chromosome");
            int positionIndex = fieldsIndexMap.get("Position");
            int refIndex = fieldsIndexMap.get("Ref");
            int altIndex = fieldsIndexMap.get("Alt");
            int codingEffectIndex = fieldsIndexMap.get("CodingEffect");
            int hgvsProteinIndex = fieldsIndexMap.get("HgvsProtein");
            int clinvarSignificanceIndex = fieldsIndexMap.get("ClinvarSignificance");
            int clinvarSigInfoIndex = fieldsIndexMap.get("ClinvarSigInfo");

            String line;

            while ((line = file.readLine()) != null)
            {
                if (line.isEmpty())
                    break;

                ++lineIndex;

                // parse CSV data. -1 as 2nd param makes sure we include trailing empty fields
                String[] items = line.split(",", -1);

                if (items.length < BACHELOR_FILTER_CSV_FIELD_COUNT)
                {
                    LOGGER.error("Invalid item count({}), lineIndex({})", items.length, lineIndex);
                    continue;
                }

                // Gene,TranscriptId,Chromosome,Position,Ref,Alt,CodingEffect,AllEffects,HgvsProtein,HgvsCoding,DBSnpId,ClinvarSignificance,ClinvarSigInfo
                VariantFilter filter = new VariantFilter(
                        items[geneIndex],
                        items[transcriptIndex],
                        items[chromosomeIndex],
                        Long.parseLong(items[positionIndex]),
                        items[refIndex],
                        items[altIndex],
                        CodingEffect.valueOf(items[codingEffectIndex]),
                        items[hgvsProteinIndex],
                        items[clinvarSignificanceIndex],
                        items[clinvarSigInfoIndex],
                        -1);

                filters.add(filter);
            }

        }
        catch (IOException e)
        {
            LOGGER.error("failed to read bachelor input CSV file({}) index({}): {}", filterFile, lineIndex, e.toString());
        }

        LOGGER.info("loaded {} ClinVar filter records from {}", filters.size(), filterFile);

        return filters;
    }

    private void createFilterFile(final Program program, final String outputDir)
    {
        if (!Files.exists(Paths.get(mInputFilterFile)))
        {
            LOGGER.error("Failed to load filter input file: {}", mInputFilterFile);
            return;
        }

        if(!initialiseFilterWriter(outputDir))
            return;

        ProgramPanel programConfig = program.getPanel().get(0);

        List<GeneIdentifier> genesList = programConfig.getGene();

        // take up a collection of the effects to search for
        mRequiredEffects = programConfig.getSnpEffect().stream().map(SnpEffect::value).collect(Collectors.toList());
        mPanelTranscripts = genesList.stream().map(GeneIdentifier::getEnsembl).collect(Collectors.toList());

        File vcfFile = new File(mInputFilterFile);

        VCFFileReader vcfReader = new VCFFileReader(vcfFile, false);

        for (VariantContext variant : vcfReader)
        {
            processVariant(variant);
        }

        closeBufferedWriter(mFilterWriter);
    }

    // Clinvar annotations
    private static final String CLINVAR_SIGNIFICANCE = "CLNSIG";
    private static final String CLINVAR_SIG_INFO = "CLNSIGCONF";
    private static final String CLINVAR_DISEASE_NAME = "CLNDN";
    private static final String CLINVAR_MC = "MC";
    private static final String CLINVAR_RS_DB_SNP_ID = "RS";
    private static final String CLINVAR_PATHOGENIC = "Pathogenic";
    private static final String CLINVAR_LIKELY_PATHOGENIC = "Likely_pathogenic";
    private static final String CLINVAR_BENIGN = "Benign";
    private static final String CLINVAR_LIKELY_BENIGN = "Likely_benign";
    private static final String CLINVAR_CONFLICTING = "Conflicting";

    private static boolean isPathogenic(final String clinvarSignificance)
    {
        return clinvarSignificance.contains(CLINVAR_PATHOGENIC) || clinvarSignificance.contains(CLINVAR_LIKELY_PATHOGENIC);
    }

    static boolean isBenign(final String clinvarSignificance)
    {
        return clinvarSignificance.contains(CLINVAR_BENIGN) || clinvarSignificance.contains(CLINVAR_LIKELY_BENIGN);
    }

    private static boolean isConflicting(final String clinvarSignificance)
    {
        return clinvarSignificance.contains(CLINVAR_CONFLICTING);
    }

    static String stripTranscriptVersion(final String transcript)
    {
        // necessary since SnpEff adds the version id in HG38
        if(transcript.contains("."))
            return transcript.substring(0, transcript.indexOf('.'));
        else
            return transcript;
    }

    private void processVariant(final VariantContext variant)
    {
        // LOGGER.debug("read var({}) chr({}))", variant.getID(), variant.getContig());

        // first check against the genes list
        final List<SnpEffAnnotation> variantAnnotations = SnpEffAnnotationFactory.fromContext(variant);

        for (SnpEffAnnotation snpEff : variantAnnotations)
        {
            if (!snpEff.isTranscriptFeature())
                continue;

            if (!mPanelTranscripts.contains(stripTranscriptVersion(snpEff.transcript())))
                continue;

            final String clinvarSignificance = stripArrayChars(variant.getCommonInfo().getAttributeAsString(CLINVAR_SIGNIFICANCE, ""));
            final String clinvarSigInfo = stripArrayChars(variant.getCommonInfo().getAttributeAsString(CLINVAR_SIG_INFO, ""));
            final String gene = snpEff.gene();
            final CodingEffect codingEffect = CodingEffect.effect(gene, snpEff.consequences());

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
        final String transcriptId = stripTranscriptVersion(snpEff.transcript());
        final String chromosome = variant.getContig();
        long position = variant.getStart();
        String ref = variant.getReference().getBaseString();
        String alt = variant.getAlleles().get(1).getBaseString();

        ref = ref.replaceAll("\\*", "");
        alt = alt.replaceAll("\\*", "");

        final String effects = snpEff.effects();
        final String clinvarDisease = variant.getCommonInfo().getAttributeAsString(CLINVAR_DISEASE_NAME, "").replaceAll(",", ";");
        final String clinvarEffects = variant.getCommonInfo().getAttributeAsString(CLINVAR_MC, "").replaceAll(",", ";");

        final String hgvsp = snpEff.hgvsProtein();
        final String hgvsc = snpEff.hgvsCoding();

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

            mFilterWriter.write(String.format(",%s,%s", hgvsp, hgvsc));

            mFilterWriter.write(String.format(",%s,%s,%s,%s",
                    clinvarSignificance, clinvarSigInfo, clinvarDisease, clinvarEffects));

            mFilterWriter.newLine();
        }
        catch(IOException e)
        {
            LOGGER.error("Error writing filter output: {}", e.toString());
        }
    }

    private boolean initialiseFilterWriter(final String outputDir)
    {
        try
        {
            String filterFileName = outputDir;

            if (!filterFileName.endsWith(File.separator))
                filterFileName += File.separator;

            filterFileName += "bachelor_clinvar_filters.csv";

            mFilterWriter = createBufferedWriter(filterFileName, false);

            mFilterWriter.write("Gene,TranscriptId,Chromosome,Position,Ref,Alt,CodingEffect,AllEffects");
            mFilterWriter.write(",HgvsProtein,HgvsCoding,ClinvarSignificance,ClinvarSigInfo,ClinvarDisease,ClinvarEffects");
            mFilterWriter.newLine();
        }
        catch (IOException e)
        {
            LOGGER.error("Failed to create output filter file: {}", e.toString());
            return false;
        }

        return true;
    }
}
