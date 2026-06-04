package com.hartwig.hmftools.sage.vis;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.util.Map.entry;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.vis.BaseSeqViewModel.fromStr;
import static com.hartwig.hmftools.common.vis.ColorUtil.DARK_BLUE;
import static com.hartwig.hmftools.sage.vis.SageVisConstants.MAX_READ_UPPER_LIMIT;
import static com.hartwig.hmftools.sage.vis.VisFileBuilder.contextViewModelFromVariant;

import java.awt.Color;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.StringJoiner;
import java.util.concurrent.ThreadLocalRandom;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.vis.BaseSeqViewModel;
import com.hartwig.hmftools.common.vis.SvgRender;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.common.ReadContextMatch;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.quality.QualityScores;
import com.hartwig.hmftools.sage.sync.FragmentData;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

public class VariantVis
{
    protected final VisConfig mConfig;
    protected final String SampleId;
    protected final VariantTier Tier;
    protected final VariantReadContext ReadContext;
    protected final BaseRegion ViewRegion;
    protected final Map<Integer, List<SvgRender.BoxBorder>> ContextBorders;
    protected final BaseSeqViewModel RefViewModel;
    protected final BaseSeqViewModel ContextViewModel;
    protected final String VariantKey;
    protected final String IndexedBasesKey;
    protected final SimpleVariant VariantInfo;

    private int mReadCount;
    private final EnumMap<ReadContextMatch, List<ReadEvidenceRecord>> mReadEvidenceRecordsByType;
    private final EnumMap<ReadContextMatch, Integer> mReadCountByType;

    protected static final Map<String, SortedSet<String>> VARIANT_INDEXED_BASES_MAP = Maps.newConcurrentMap();

    public VariantVis(
            final SageConfig config, final String sample, final SimpleVariant variant, final VariantReadContext readContext,
            final VariantTier variantTier)
    {
        mConfig = config.Visualiser;
        SampleId = sample;
        VariantInfo = variant;
        Tier = variantTier;
        mReadEvidenceRecordsByType = Maps.newEnumMap(ReadContextMatch.class);
        ReadContext = readContext;
        VariantKey = VariantInfo.chromosome() + "_" + VariantInfo.position() + "_" + VariantInfo.ref() + "_" + VariantInfo.alt();

        int coreStart = ReadContext.CorePositionStart;
        int coreEnd = ReadContext.CorePositionEnd;
        int flankStart = coreStart - ReadContext.leftFlankLength();
        int flankEnd = coreEnd + ReadContext.rightFlankLength();
        ViewRegion = new BaseRegion(coreStart - SageVisConstants.READ_EXTEND_LENGTH, coreEnd + SageVisConstants.READ_EXTEND_LENGTH);

        ContextBorders = Map.ofEntries(
                entry(VariantInfo.Position,
                        Lists.newArrayList(
                                new SvgRender.BoxBorder(SvgRender.BorderLocation.LEFT, Color.BLACK),
                                new SvgRender.BoxBorder(SvgRender.BorderLocation.RIGHT, Color.BLACK))),
                entry(coreStart,
                        Lists.newArrayList(
                                new SvgRender.BoxBorder(SvgRender.BorderLocation.LEFT, Color.BLUE))),
                entry(coreEnd,
                        Lists.newArrayList(
                                new SvgRender.BoxBorder(SvgRender.BorderLocation.RIGHT, Color.BLUE))),
                entry(flankStart,
                        Lists.newArrayList(
                                new SvgRender.BoxBorder(SvgRender.BorderLocation.LEFT, DARK_BLUE))),
                entry(flankEnd,
                        Lists.newArrayList(
                                new SvgRender.BoxBorder(SvgRender.BorderLocation.RIGHT, DARK_BLUE)))
        );

        RefGenomeCoordinates refGenCoords = config.RefGenVersion == V37 ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;
        int chromosomeLength = refGenCoords.length(VariantInfo.chromosome());
        int refPosStart = max(1, ViewRegion.start());
        int refPosEnd = min(chromosomeLength, ViewRegion.end());

        RefGenomeSource refGenome = loadRefGenome(config.RefGenomeFile);
        String refBases = refGenome.getBaseString(VariantInfo.chromosome(), refPosStart, refPosEnd);
        RefViewModel = fromStr(refBases, refPosStart);

        ContextViewModel = contextViewModelFromVariant(ReadContext, VariantInfo.ref(), VariantInfo.alt());

        StringJoiner indexedBasesKeyBuilder = new StringJoiner("_");
        indexedBasesKeyBuilder.add(String.valueOf(ReadContext.VarIndex));
        indexedBasesKeyBuilder.add(String.valueOf(ReadContext.CoreIndexStart));
        indexedBasesKeyBuilder.add(String.valueOf(ReadContext.CoreIndexEnd));
        indexedBasesKeyBuilder.add(String.valueOf(ReadContext.leftFlankLength()));
        indexedBasesKeyBuilder.add(new String(ReadContext.ReadBases));
        IndexedBasesKey = indexedBasesKeyBuilder.toString();

        mReadCountByType = Maps.newEnumMap(ReadContextMatch.class);
        mReadCount = 0;
        VARIANT_INDEXED_BASES_MAP.computeIfAbsent(VariantKey, k -> Sets.newTreeSet()).add(IndexedBasesKey);
    }

    public int readCount() { return mReadCount; }
    public EnumMap<ReadContextMatch, List<ReadEvidenceRecord>> readEvidenceRecordsByType() { return mReadEvidenceRecordsByType; }
    public EnumMap<ReadContextMatch, Integer> readCountByType() { return mReadCountByType; }

    public void addEvidence(
            final SAMRecord read, @Nullable final FragmentData fragment, final ReadContextMatch matchType,
            @Nullable final QualityScores modifiedQualities)
    {
        ++mReadCount;
        mReadCountByType.put(matchType, mReadCountByType.getOrDefault(matchType, 0) + 1);

        List<ReadEvidenceRecord> records = mReadEvidenceRecordsByType.computeIfAbsent(matchType, k -> Lists.newArrayList());
        if(fragment == null)
        {
            records.add(new ReadEvidenceRecord(read, null, matchType, modifiedQualities, VariantInfo.Position));
            return;
        }

        BaseRegion firstUnclippedRegion = new BaseRegion(fragment.First.getUnclippedStart(), fragment.First.getUnclippedEnd());
        boolean firstIsVisible = ViewRegion.overlaps(firstUnclippedRegion);

        BaseRegion secondUnclippedRegion = new BaseRegion(fragment.Second.getUnclippedStart(), fragment.Second.getUnclippedEnd());
        boolean secondIsVisible = ViewRegion.overlaps(secondUnclippedRegion);
        if(firstIsVisible && !secondIsVisible)
        {
            records.add(new ReadEvidenceRecord(fragment.First, null, matchType, modifiedQualities, VariantInfo.Position));
            return;
        }

        if(!firstIsVisible && secondIsVisible)
        {
            records.add(new ReadEvidenceRecord(fragment.Second, null, matchType, modifiedQualities, VariantInfo.Position));
            return;
        }

        records.add(new ReadEvidenceRecord(read, fragment, matchType, modifiedQualities, VariantInfo.Position));
    }

    protected void downsampleReadEvidenceRecords()
    {
        for(Map.Entry<ReadContextMatch, List<ReadEvidenceRecord>> entry : mReadEvidenceRecordsByType.entrySet())
        {
            ReadContextMatch matchType = entry.getKey();
            List<ReadEvidenceRecord> records = entry.getValue();

            int maxReads = SageVisConstants.MAX_READS_PER_TYPE.get(matchType);

            if(mConfig.MaxSupportReads < 0)
            {
                maxReads = MAX_READ_UPPER_LIMIT;
            }
            else if(mConfig.MaxSupportReads > 0)
            {
                maxReads = min(mConfig.MaxSupportReads, MAX_READ_UPPER_LIMIT);
            }

            List<ReadEvidenceRecord> newRecords = Lists.newArrayList();
            while(!records.isEmpty() && newRecords.size() < maxReads)
            {
                if(records.size() == 1)
                {
                    newRecords.add(records.get(0));
                    records.remove(0);
                    continue;
                }

                int i = ThreadLocalRandom.current().nextInt(records.size());
                if(i < records.size() - 1)
                {
                    ReadEvidenceRecord tmp = records.get(i);
                    records.set(i, records.get(records.size() - 1));
                    records.set(records.size() - 1, tmp);
                }

                newRecords.add(records.get(records.size() - 1));
                records.remove(records.size() - 1);
            }

            mReadEvidenceRecordsByType.put(matchType, newRecords);
        }
    }
}
