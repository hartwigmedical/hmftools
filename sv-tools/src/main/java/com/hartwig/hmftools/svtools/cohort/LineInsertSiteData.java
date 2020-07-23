package com.hartwig.hmftools.svtools.cohort;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionsOverlap;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.utils.sv.SvRegion;

import org.apache.commons.compress.utils.Lists;

public class LineInsertSiteData
{
    public final String Program;
    public final String SampleId;

    public final SvRegion InsertSite;
    public final String InsertType;
    public final SvRegion SourceSite;
    public final boolean HasInversion;

    // Linx only
    public final String ChainDesc;
    public final int ClusterId;
    public final int ChainId;

    public static final String INSERT_TYPE_SOLO_L1 = "SOLO_L1";
    public static final String INSERT_TYPE_PARTNERED = "PARTNERED";
    public static final String INSERT_TYPE_TRANSDUCTION = "TRANSDUCTION";

    public static final String PROGRAM_LINX = "LINX";
    public static final String PROGRAM_PCAWG = "PCAWG";

    public LineInsertSiteData(final String program, final String sampleId, final SvRegion insertSite, final String insertType,
            final SvRegion sourceSite, final boolean hasInversion, final int clusterId, final int chainId, final String chainDesc)
    {
        Program = program;
        SampleId = sampleId;
        InsertSite = insertSite;
        InsertType = insertType;
        SourceSite = sourceSite;
        HasInversion = hasInversion;
        ChainDesc = chainDesc;
        ClusterId = clusterId;
        ChainId = chainId;
    }

    private static final int MAX_BASE_DIFF = 50;

    public boolean matches(final LineInsertSiteData other)
    {
        if(!other.InsertSite.Chromosome.equals(InsertSite.Chromosome))
            return false;

        if(other.InsertSite.isValid() && InsertSite.isValid())
        {
            return positionsOverlap(
                    other.InsertSite.start(), other.InsertSite.end(),
                    InsertSite.start() - MAX_BASE_DIFF, InsertSite.end() + MAX_BASE_DIFF);
        }

        // either of the positions may be set
        final int[] pos1 = {
                InsertSite.start() > 0 ? InsertSite.start() : InsertSite.end(),
                InsertSite.end() > 0 ? InsertSite.end() : InsertSite.start() };

        final int[] pos2 = {
                other.InsertSite.start() > 0 ? other.InsertSite.start() : other.InsertSite.end(),
                other.InsertSite.end() > 0 ? other.InsertSite.end() : other.InsertSite.start() };

        return positionsOverlap(pos1[SE_START], pos1[SE_END], pos2[SE_START] - MAX_BASE_DIFF, pos2[SE_END] + MAX_BASE_DIFF);
    }

    public static LineInsertSiteData fromExternal(final Map<String, Integer> fieldsIndexMap, final String data)
    {
        final String[] items = data.split(",", -1);

        final String sampleId = items[fieldsIndexMap.get("SampleId")];

        // SampleId,Chromosome,Pos1,Pos2,InsertType,SourceSite,TransductionSite

        final String insertTypeRaw = items[fieldsIndexMap.get("InsertType")];
        final String sourceSiteStr = items[fieldsIndexMap.get("SourceSite")];

        final String[] sourcePosItems = sourceSiteStr.split("_");

        // 22_29059272_29065304
        final SvRegion sourceRegion = sourcePosItems.length == 3 ?
                new SvRegion(sourcePosItems[0], Integer.parseInt(sourcePosItems[1]), Integer.parseInt(sourcePosItems[2])) : null;

        final String insertChromosome = items[fieldsIndexMap.get("InsertChr")];

        final int[] insertPositions = new int[] { Integer.parseInt(items[fieldsIndexMap.get("InsertPos1")]), -1 };
        final String insertPos2Str = items[fieldsIndexMap.get("InsertPos2")];

        int lowerInsertIndex = 0;
        if(!insertPos2Str.equals("UNK"))
        {
            insertPositions[1] = Integer.parseInt(insertPos2Str);
            lowerInsertIndex = insertPositions[0] <= insertPositions[1] ? 0 : 1;
        }

        if(insertChromosome.isEmpty() || insertChromosome.equals("0"))
            return null;

        final SvRegion insertRegion = new SvRegion(
                insertChromosome, insertPositions[lowerInsertIndex], insertPositions[switchIndex(lowerInsertIndex)]);

        final String insertType = insertTypeRaw.equals("Solo-L1") ? INSERT_TYPE_SOLO_L1 :
                (insertTypeRaw.equals("Orphan") ? INSERT_TYPE_TRANSDUCTION : INSERT_TYPE_PARTNERED);

        final boolean hasInversion = items[fieldsIndexMap.get("InsertStructure")].equals("INV");

        return new LineInsertSiteData(
                PROGRAM_PCAWG, sampleId, insertRegion, insertType, sourceRegion, hasInversion, -1, -1, "");
    }

    public static LineInsertSiteData fromLinx(final Map<String, Integer> fieldsIndexMap, final String data)
    {
        final String[] items = data.split(",", -1);

        final String sampleId = items[fieldsIndexMap.get("SampleId")];
        final int clusterId = Integer.parseInt(items[fieldsIndexMap.get("ClusterId")]);
        final int chainId = Integer.parseInt(items[fieldsIndexMap.get("ChainId")]);
        final String chainDesc = items[fieldsIndexMap.get("ChainDesc")];

        final String insertChromosome = items[fieldsIndexMap.get("InsertChr")];

        final int[] insertPositions =
                new int[] { Integer.parseInt(items[fieldsIndexMap.get("InsertPosStart")]),
                        Integer.parseInt(items[fieldsIndexMap.get("InsertPosEnd")]) };

        if(insertChromosome.isEmpty() || insertChromosome.equals("0"))
            return null;

        final String sourceChromosome = items[fieldsIndexMap.get("SourceChr")];

        final int[] sourcePositions =
                new int[] { Integer.parseInt(items[fieldsIndexMap.get("SourcePosStart")]),
                        Integer.parseInt(items[fieldsIndexMap.get("SourcePosEnd")]) };

        final int[] invPositions =
                new int[] { Integer.parseInt(items[fieldsIndexMap.get("SourceInvPosStart")]),
                        Integer.parseInt(items[fieldsIndexMap.get("SourceInvPosEnd")]) };


        final SvRegion insertRegion = new SvRegion(insertChromosome, insertPositions);
        final SvRegion sourceRegion = !sourceChromosome.isEmpty() ? new SvRegion(sourceChromosome, sourcePositions) : null;

        final String insertType = sourceRegion == null ? INSERT_TYPE_SOLO_L1 : INSERT_TYPE_TRANSDUCTION;

        final boolean hasInversion = invPositions[SE_START] > 0;

        return new LineInsertSiteData(
                PROGRAM_LINX, sampleId, insertRegion, insertType, sourceRegion, hasInversion, clusterId, chainId, chainDesc);

    }
}