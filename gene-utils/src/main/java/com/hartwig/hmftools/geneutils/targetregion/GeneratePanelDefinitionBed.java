package com.hartwig.hmftools.geneutils.targetregion;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.ENHANCER;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_5P;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.SpecificRegions.parseStandardFormat;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;
import static com.hartwig.hmftools.geneutils.targetregion.GenerateTargetRegionsBed.OUTPUT_FILE;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptCodingType;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;

public class GeneratePanelDefinitionBed
{
    private final RefGenomeVersion mRefGenVersion;
    private final EnsemblDataCache mEnsemblDataCache;

    private final String mOutputFile;

    private final List<RegionDefinition> mRegions;
    private final List<ChrBaseRegion> mCombinedRegions;

    private static final String PANEL_DEFINITION_FILE = "panel_definition_file";

    private static final int LOCATION_PROBE_BUFFER = 120;
    private static final int PROMOTER_PROBE_BUFFER = 1000;

    public GeneratePanelDefinitionBed(final ConfigBuilder configBuilder)
    {
        mRefGenVersion = RefGenomeVersion.from(configBuilder);
        mEnsemblDataCache = new EnsemblDataCache(configBuilder);
        mEnsemblDataCache.setRequiredData(true, false, false, true);
        mEnsemblDataCache.load(false);

        mEnsemblDataCache.createGeneNameIdMap();

        String panelDefintion = configBuilder.getValue(PANEL_DEFINITION_FILE);
        mRegions = loadDefinitionRegions(panelDefintion);

        mOutputFile = configBuilder.getValue(OUTPUT_FILE);
        mCombinedRegions = Lists.newArrayList();
    }

    private enum RegionType
    {
        LOCATION,
        GENE;
    }

    private class RegionDefinition
    {
        public final RegionType Type;
        public final String Label;
        public final String GeneName;
        public final List<TranscriptCodingType> CodingTypes;
        public final List<Integer> Exons;

        public final ChrBaseRegion SpecificRegion;

        public RegionDefinition(
                final String label, final RegionType type,
                final String geneName, final List<TranscriptCodingType> codingTypes, final List<Integer> exons)
        {
            Label = label;
            Type = type;
            GeneName = geneName;
            CodingTypes = codingTypes;
            Exons = exons;
            SpecificRegion = null;
        }

        public RegionDefinition(
                final String label, final RegionType type, ChrBaseRegion specificRegion)
        {
            Label = label;
            Type = type;
            SpecificRegion = specificRegion;
            GeneName = "";
            CodingTypes = Collections.emptyList();
            Exons = Collections.emptyList();
        }
    }

    private List<RegionDefinition> loadDefinitionRegions(final String filename)
    {
        if(filename == null)
            return Collections.emptyList();

        List<RegionDefinition> regions = Lists.newArrayList();
        try
        {
            final List<String> fileContents = Files.readAllLines(new File(filename).toPath());

            fileContents.remove(0);

            for(final String line : fileContents)
            {
                String[] values = line.split(TSV_DELIM, -1);
                try
                {
                    RegionType regionType = RegionType.valueOf(values[0]);

                    if(regionType == RegionType.LOCATION)
                    {
                        String[] locationInfo = values[3].split(":");
                        String chromosome = locationInfo[0];
                        int position = Integer.parseInt(locationInfo[1]);

                        ChrBaseRegion specificRegion = new ChrBaseRegion(
                                chromosome, position - LOCATION_PROBE_BUFFER, position + LOCATION_PROBE_BUFFER);

                        String label = values[2];

                        regions.add(new RegionDefinition(label, regionType, specificRegion));
                    }
                    else
                    {
                        String gene = values[1];
                        String label = values[2];
                        String regionStr = values[3];

                        if(regionStr.startsWith("region"))
                        {
                            String regionInfo = regionStr.split("=", 2)[1];

                            ChrBaseRegion specificRegion = parseStandardFormat(regionInfo);
                            regions.add(new RegionDefinition(gene, RegionType.LOCATION, specificRegion));
                            continue;
                        }

                        List<TranscriptCodingType> codingTypes = Lists.newArrayList();
                        List<Integer> exons = Lists.newArrayList();

                        if(regionStr.equals("CDS"))
                        {
                            codingTypes.add(TranscriptCodingType.CODING);
                        }
                        else if(regionStr.equals("5'UTR"))
                        {
                            codingTypes.add(UTR_5P);
                        }
                        else if(regionStr.equals("3'UTR"))
                        {
                            codingTypes.add(TranscriptCodingType.UTR_3P);
                        }
                        else if(regionStr.equals("Promoter"))
                        {
                            codingTypes.add(TranscriptCodingType.ENHANCER);
                        }
                        else if(regionStr.contains("ex"))
                        {
                            parseExonData(exons, regionStr);
                        }
                        else
                        {
                            GU_LOGGER.warn("unhandled region type({}) from line({})", regionStr, line);
                            continue;
                        }

                        regions.add(new RegionDefinition(label, regionType, gene, codingTypes, exons));
                    }
                }
                catch(Exception e)
                {
                    GU_LOGGER.error("error parsing region definition data({}): {}", line, e.toString());
                }
            }

            GU_LOGGER.info("loaded {} panel definition regions from file: {}", regions.size(), filename);
            return regions;
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to loadpanel definition regions from file({}): {}", filename, e.toString());
            return Collections.emptyList();
        }
    }

    private static void parseExonData(final List<Integer> exons, final String exonDataStr)
    {
        // like ex13-17,24
        String[] exonStrList = exonDataStr.substring(2).split(",");

        for(String exonStr : exonStrList)
        {
            if(exonStr.contains("-"))
            {
                String[] exonStartEnd = exonStr.split("-", 2);

                int exonStart = Integer.parseInt(exonStartEnd[0]);
                int exonEnd = Integer.parseInt(exonStartEnd[1]);

                for(int j = exonStart; j <= exonEnd; ++j)
                {
                    exons.add(j);
                }
            }
            else
            {
                exons.add(Integer.parseInt(exonStr));
            }
        }
    }

    public void run()
    {
        if(mRegions.isEmpty() || mOutputFile == null)
        {
            GU_LOGGER.error("invalid config or input file, exiting");
            System.exit(1);
        }

        GU_LOGGER.info("generating BED file from {} panel definition regions to {}", mRegions.size(), mOutputFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputFile, false);

            // no header
            List<LabelledRegion> labelledRegions = Lists.newArrayList();

            for(RegionDefinition region : mRegions)
            {
                if(region.Type == RegionType.LOCATION)
                {
                    labelledRegions.add(buildRegion(
                            region.SpecificRegion.Chromosome, region.SpecificRegion.start(), region.SpecificRegion.end(),
                            region.Label, false));
                }
                else
                {
                    labelledRegions.addAll(buildGeneDataRegions(writer, region));
                }
            }

            Collections.sort(labelledRegions);

            for(LabelledRegion region : labelledRegions)
            {
                writer.write(format("%s\t%d\t%d\t%s", region.chromosome(), region.start(), region.end(), region.Label));
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to write to line ref-genome bases: {}", e.toString());
        }

        GU_LOGGER.info("panel definition BED file generation complete");
    }

    private class LabelledRegion extends ChrBaseRegion
    {
        public final String Label;

        public LabelledRegion(final String chromosome, final int posStart, final int posEnd, final String label)
        {
            super(chromosome, posStart, posEnd);
            Label = label;
        }
    }

    private LabelledRegion buildRegion(
            final String chromosome, final int posStart, final int posEnd, final String label, final boolean isCoding)
    {
        return new LabelledRegion(chromosome, posStart, posEnd, isCoding ? label + "_CODING" : label);
    }

    private List<LabelledRegion> buildGeneDataRegions(final BufferedWriter writer, final RegionDefinition region) throws IOException
    {
        GeneData geneData = mEnsemblDataCache.getGeneDataByName(region.GeneName);

        if(geneData == null)
        {
            GU_LOGGER.error("unknown gene({}) label({})", region.GeneName, region.Label);
            Collections.emptyList();
        }

        List<TranscriptData> transcripts = mEnsemblDataCache.getTranscripts(geneData.GeneId);

        List<LabelledRegion> regions = Lists.newArrayList();

        for(TranscriptData transData : transcripts)
        {
            if(transData.CodingStart == null || !transData.IsCanonical)
                continue;

            if(region.CodingTypes.contains(ENHANCER))
            {
                int promoterStart, promoterEnd;

                if(transData.Strand == POS_ORIENT)
                {
                    promoterStart = transData.TransStart - PROMOTER_PROBE_BUFFER;
                    promoterEnd = transData.TransStart;
                }
                else
                {
                    promoterStart = transData.TransEnd;
                    promoterEnd = transData.TransEnd + PROMOTER_PROBE_BUFFER;
                }

                regions.add(buildRegion(
                        geneData.Chromosome, promoterStart, promoterEnd, format("%s_PROMOTER", geneData.GeneName), false));

                if(region.CodingTypes.size() == 1)
                    return regions;
            }

            for(ExonData exon : transData.exons())
            {
                boolean hasUtr = positionWithin(transData.CodingStart, exon.Start, exon.End)
                        || positionWithin(transData.CodingEnd, exon.Start, exon.End)
                        || exon.End < transData.CodingStart || exon.Start > transData.CodingEnd;

                boolean hasCoding = positionsOverlap(exon.Start, exon.End, transData.CodingStart, transData.CodingEnd);

                if(hasUtr)
                {
                    TranscriptCodingType utrType;

                    if(transData.Strand == POS_ORIENT)
                        utrType = exon.Start < transData.CodingStart ? UTR_5P : UTR_3P;
                    else
                        utrType = exon.End > transData.CodingEnd ? UTR_5P : UTR_3P;

                    if(region.CodingTypes.contains(utrType))
                    {
                        int utrStart, utrEnd;

                        if(exon.Start < transData.CodingStart)
                        {
                            utrStart = exon.Start;
                            utrEnd = min(exon.End, transData.CodingStart);
                        }
                        else
                        {
                            utrStart = max(exon.Start, transData.CodingEnd);
                            utrEnd = exon.End;
                        }

                        regions.add(buildRegion(
                                geneData.Chromosome, utrStart, utrEnd, format("%s_%s", geneData.GeneName, utrType), false));
                    }
                }

                if(hasCoding)
                {
                    if(region.CodingTypes.contains(TranscriptCodingType.CODING) || region.Exons.contains(exon.Rank))
                    {
                        int codingStart = max(transData.CodingStart, exon.Start);
                        int codingEnd = min(transData.CodingEnd, exon.End);

                        regions.add(buildRegion(
                                geneData.Chromosome, codingStart, codingEnd, format("%s_exon%d", geneData.GeneName, exon.Rank), true));
                    }
                }
            }
        }

        return regions;
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addPath(PANEL_DEFINITION_FILE, true,"Additional regions beyond panel definition BED");
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output BED filename");

        addEnsemblDir(configBuilder, true);
        addRefGenomeVersion(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        GeneratePanelDefinitionBed panelDefinitionBed = new GeneratePanelDefinitionBed(configBuilder);
        panelDefinitionBed.run();
    }
}
