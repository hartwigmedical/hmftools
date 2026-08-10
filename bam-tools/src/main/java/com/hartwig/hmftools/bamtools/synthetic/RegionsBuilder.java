package com.hartwig.hmftools.bamtools.synthetic;

import static java.lang.String.format;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.bamtools.synthetic.RegionType.CN_BACKBONE;
import static com.hartwig.hmftools.bamtools.synthetic.RegionType.DISRUPTION;
import static com.hartwig.hmftools.bamtools.synthetic.RegionType.FUSION;
import static com.hartwig.hmftools.bamtools.synthetic.RegionType.GERMLINE_MUTATION;
import static com.hartwig.hmftools.bamtools.synthetic.RegionType.GERMLINE_SV;
import static com.hartwig.hmftools.bamtools.synthetic.RegionType.PANEL;
import static com.hartwig.hmftools.bamtools.synthetic.RegionType.SOMATIC_MUTATION;
import static com.hartwig.hmftools.bamtools.synthetic.RegionType.SOMATIC_SV;
import static com.hartwig.hmftools.bamtools.synthetic.RegionType.typesStr;
import static com.hartwig.hmftools.bamtools.synthetic.SyntheticConstants.SMALL_VARIANT_BUFFER;
import static com.hartwig.hmftools.bamtools.synthetic.SyntheticConstants.SV_DISCORDANT_BUFFER;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.region.ChrBaseRegion.loadChrBaseRegionList;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INF;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REGION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_REGION_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.variant.SmallVariant;
import com.hartwig.hmftools.common.variant.SmallVariantFactory;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;

public class RegionsBuilder
{
    private final SyntheticConfig mConfig;

    private final List<RegionData> mRegions;
    private List<PurpleCopyNumber> mCopyNumbers;

    public RegionsBuilder(final SyntheticConfig config)
    {
        mConfig = config;
        mRegions = Lists.newArrayList();
        mCopyNumbers = Lists.newArrayList();
    }

    public List<RegionData> regions() { return mRegions; }

    public void buildRegions()
    {
        if(mConfig.RegionsFilename != null)
        {
            loadPanelRegions();
        }

        if(mConfig.PurpleDir != null)
        {
            if(!loadPurpleRegions())
                System.exit(1);
        }

        if(mConfig.LinxDir != null)
        {
            if(!loadLinxRegions())
                System.exit(1);
        }

        addCopyNumberBackbone();

        mergeRegions();

        setRegionCopyNumbers();

        logRegionStats();

        if(mConfig.WriteRegionData)
            writeRegionData();
    }

    private boolean loadPanelRegions()
    {
        List<ChrBaseRegion> regions = loadChrBaseRegionList(mConfig.RegionsFilename);

        for(ChrBaseRegion region : regions)
        {
            RegionData regionData = new RegionData(region.chromosome(), region.start(), region.end(), PANEL);
            mRegions.add(regionData);
        }

        return true;
    }

    private void addCopyNumberBackbone()
    {
        if(!mConfig.isWgsMode() || mConfig.CopyNumberBackboneDistance <= 0)
            return;

        RefGenomeCoordinates refGenomeCoordinates = RefGenomeCoordinates.refGenomeCoordinates(mConfig.RefGenVersion);

        for(Map.Entry<Chromosome,Integer> entry : refGenomeCoordinates.Lengths.entrySet())
        {
            String chrStr = mConfig.RefGenVersion.versionedChromosome(entry.getKey().toString());
            int length = entry.getValue();

            for(int position = mConfig.CopyNumberBackboneDistance; position < length; position += mConfig.CopyNumberBackboneDistance)
            {
                RegionData regionData = new RegionData(chrStr, position, position, CN_BACKBONE);
                mRegions.add(regionData);
            }
        }
    }

    private boolean loadPurpleRegions()
    {
        try
        {
            String copyNumberFile = PurpleCopyNumberFile.generateFilenameForReading(mConfig.PurpleDir, mConfig.SampleId);
            mCopyNumbers.addAll(PurpleCopyNumberFile.read(copyNumberFile));

            addSmallVariants(true);

            if(mConfig.ReferenceId != null)
                addSmallVariants(false);

            addSomaticStructuralVariants();
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to to load Purple data: {}", e.toString());
            return false;
        }

        return true;
    }

    private void addSmallVariants(boolean isSomatic) throws IOException
    {
        String vcfFile = isSomatic ?
                PurpleCommon.purpleSomaticVcfFile(mConfig.PurpleDir, mConfig.SampleId)
                : PurpleCommon.purpleGermlineVcfFile(mConfig.PurpleDir, mConfig.SampleId);

        List<SmallVariant> variants = SmallVariantFactory.passOnlyInstance().fromVCFFile(
                mConfig.SampleId, mConfig.ReferenceId, null, vcfFile);

        int count = 0;

        for(SmallVariant variant : variants)
        {
            if(!variant.reported())
                continue;

            RegionData regionData = new RegionData(
                    variant.chromosome(),
                    variant.position() - SMALL_VARIANT_BUFFER,
                    variant.position() + SMALL_VARIANT_BUFFER,
                    isSomatic ? SOMATIC_MUTATION : GERMLINE_MUTATION);

            mRegions.add(regionData);

            ++count;

            if(count > mConfig.MaxSmallVariants)
            {
                BT_LOGGER.debug("reached max small variant count({})", mConfig.MaxSmallVariants);
                break;
            }
        }
    }

    private void addSomaticStructuralVariants() throws IOException
    {
        String vcfFile = PurpleCommon.purpleSomaticSvFile(mConfig.PurpleDir, mConfig.SampleId);

        List<StructuralVariant> variants = StructuralVariantFileLoader.fromFile(vcfFile, new AlwaysPassFilter());

        int count = 0;

        for(StructuralVariant sv : variants)
        {
            if(!sv.filter().equals(PASS_FILTER))
                continue;

            if(sv.type() == INF)
                continue;

            int positionBuffer = 0;

            if(sv.type() != SGL && applyDiscordantBuffer(sv.type(), sv.position(true), sv.position(false)))
                positionBuffer = SV_DISCORDANT_BUFFER;

            if(sv.type() == BND)

            for(int se = SE_START; se <= SE_END; ++se)
            {
                if(sv.type() == SGL && se == SE_END)
                    break;

                int position = sv.position(isStart(se));

                RegionData regionData = new RegionData(
                        sv.chromosome(isStart(se)), position - positionBuffer, position + positionBuffer, SOMATIC_SV);

                mRegions.add(regionData);
            }

            ++count;

            if(count > mConfig.MaxSVs)
            {
                BT_LOGGER.debug("reached max SV count({})", mConfig.MaxSmallVariants);
                break;
            }
        }
    }

    private static boolean applyDiscordantBuffer(final StructuralVariantType type, int posStart, int posEnd)
    {
        if(type == BND)
            return true;

        if(type == DEL || type == DUP || type == INV)
            return posEnd - posStart >= 2000;
        else
            return false;
    }

    private boolean loadLinxRegions()
    {
        try
        {
            String somaticSvAnnotationFile = LinxSvAnnotation.generateFilename(mConfig.LinxDir, mConfig.SampleId, false);
            String somaticBreakendFile = LinxBreakend.generateFilename(mConfig.LinxDir, mConfig.SampleId);
            String fusionFile = LinxFusion.generateFilename(mConfig.LinxDir, mConfig.SampleId);

            List<LinxFusion> fusions = LinxFusion.read(fusionFile);
            List<LinxSvAnnotation> somaticSvAnnotations = LinxSvAnnotation.read(somaticSvAnnotationFile);
            List<LinxBreakend> somaticBreakends = LinxBreakend.read(somaticBreakendFile);

            Set<LinxBreakend> processedBreakends = Sets.newHashSet();

            for(LinxFusion fusion : fusions)
            {
                if(fusion.reported())
                {
                    LinxBreakend breakendFive = somaticBreakends.stream()
                            .filter(x -> x.id() == fusion.fivePrimeBreakendId()).findFirst().orElse(null);

                    LinxBreakend breakendThree = somaticBreakends.stream()
                            .filter(x -> x.id() == fusion.threePrimeBreakendId()).findFirst().orElse(null);

                    addStructuralVariantRegion(breakendFive, somaticSvAnnotations, processedBreakends, FUSION);
                    addStructuralVariantRegion(breakendThree, somaticSvAnnotations, processedBreakends, FUSION);
                }
            }

            for(LinxBreakend breakend : somaticBreakends)
            {
                if(breakend.reportedStatus() == ReportedStatus.REPORTED)
                {
                    addStructuralVariantRegion(breakend, somaticSvAnnotations, processedBreakends, DISRUPTION);
                }
            }

            if(mConfig.ReferenceId != null)
            {
                String germlineSvAnnotationFile = LinxSvAnnotation.generateFilename(mConfig.LinxDir, mConfig.SampleId, true);
                String germlineBreakendFile = LinxBreakend.generateFilename(mConfig.LinxDir, mConfig.SampleId, true);
                // String germlineDisruptionFile = LinxGermlineDisruption.generateFilename(mConfig.LinxDir, mConfig.SampleId);

                List<LinxBreakend> germlineBreakends = LinxBreakend.read(germlineBreakendFile);
                List<LinxSvAnnotation> germlineSvAnnotations = LinxSvAnnotation.read(germlineSvAnnotationFile);

                for(LinxBreakend breakend : germlineBreakends)
                {
                    if(breakend.reportedStatus() == ReportedStatus.REPORTED)
                    {
                        addStructuralVariantRegion(breakend, germlineSvAnnotations, processedBreakends, GERMLINE_SV);
                    }
                }
            }
        }
        catch(IOException e)
        {
            BT_LOGGER.error("failed to to load Linx data: {}", e.toString());
            return false;
        }

        return true;
    }

    private void addStructuralVariantRegion(
            final LinxBreakend breakend, final List<LinxSvAnnotation> svAnnotations, final Set<LinxBreakend> processedBreakends,
            final RegionType regionType)
    {
        if(breakend == null || processedBreakends.contains(breakend))
            return;

        processedBreakends.add(breakend);

        LinxSvAnnotation svAnnotation = svAnnotations.stream().filter(x -> x.svId() == x.svId()).findFirst().orElse(null);

        if(svAnnotation == null)
            return;

        String chromosome = LinxBreakend.chromosomeFromCoords(breakend.coords());
        int position = LinxBreakend.positionFromCoords(breakend.coords());

        int positionBuffer = 0;

        if(svAnnotation.type() != SGL)
        {
            int positionStart = LinxBreakend.positionFromCoords(svAnnotation.coordsStart());
            int positionEnd = LinxBreakend.positionFromCoords(svAnnotation.coordsEnd());

            if(applyDiscordantBuffer(svAnnotation.type(), positionStart, positionEnd))
                positionBuffer = SV_DISCORDANT_BUFFER;
        }

        RegionData regionData = new RegionData(
                chromosome, position - positionBuffer, position + positionBuffer, regionType);

        mRegions.add(regionData);
    }

    private void mergeRegions()
    {
        Collections.sort(mRegions);

        int index = 0;
        while(index < mRegions.size() - 1)
        {
            RegionData region = mRegions.get(index);
            RegionData nextRegion = mRegions.get(index + 1);

            // merge if within the event range
            if(region.Chromosome.equals(nextRegion.Chromosome) && region.end() >= nextRegion.start() - mConfig.Params.ReadLength)
            {
                region.merge(nextRegion);
                mRegions.remove(index + 1);
            }
            else
            {
                ++index;
            }
        }
    }

    private void setRegionCopyNumbers()
    {
        if(mCopyNumbers.isEmpty())
            return;

        for(RegionData regionData : mRegions)
        {
            for(PurpleCopyNumber copyNumber : mCopyNumbers)
            {
                if(!copyNumber.chromosome().equals(regionData.Chromosome))
                    continue;

                if(positionsOverlap(regionData.start(), regionData.end(), copyNumber.start(), copyNumber.end()))
                {
                    regionData.setCopyNumber(copyNumber.averageTumorCopyNumber());
                }
            }
        }
    }

    private void logRegionStats()
    {
        long[] regionTypeLengths = new long[RegionType.values().length];
        int[] regionTypeCounts = new int[RegionType.values().length];
        long totalLength = 0;

        for(RegionData region : mRegions)
        {
            int length = region.baseLength();
            totalLength += length;

            for(RegionType type : region.types())
            {
                ++regionTypeCounts[type.ordinal()];
                regionTypeLengths[type.ordinal()] += length;
            }
        }

        BT_LOGGER.info("total slice length({})", totalLength);

        for(RegionType type : RegionType.values())
        {
            if(regionTypeLengths[type.ordinal()] != 0)
            {
                BT_LOGGER.info("region type({}) count({}) length({})",
                        type, regionTypeCounts[type.ordinal()], regionTypeLengths[type.ordinal()]);
            }
        }
    }

    private void writeRegionData()
    {
        try
        {
            String filename = mConfig.formFilename(".region_data.tsv");
            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(FLD_CHROMOSOME).add(FLD_REGION_START).add(FLD_REGION_END);
            sj.add("Types").add("CopyNumber");

            writer.write(sj.toString());
            writer.newLine();

            for(RegionData regionData : mRegions)
            {
                sj = new StringJoiner(TSV_DELIM);
                sj.add(regionData.Chromosome);
                sj.add(String.valueOf(regionData.start()));
                sj.add(String.valueOf(regionData.end()));

                sj.add(typesStr(regionData.types()));
                sj.add(format("%.2f", regionData.copyNumber()));

                writer.write(sj.toString());
                writer.newLine();
            }

            writer.close();

        }
        catch(IOException e)
        {
            BT_LOGGER.error(" failed to write region data: {}", e.toString());
        }
    }
}
