package com.hartwig.hmftools.isofox.neo;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_PAIR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.neo.AminoAcidConverter.reverseStrandBases;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import java.util.List;

import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.neo.NeoEpitopeFile;

import org.apache.commons.compress.utils.Lists;

public class NeoEpitopeData
{
    public final NeoEpitopeFile Source;

    public final String[] Chromosomes;
    public final int[] Positions;
    public final byte[] Orientations;

    public final List<String>[] Transcripts;

    public NeoEpitopeData(final NeoEpitopeFile source)
    {
        Source = source;
        Chromosomes = new String[FS_PAIR];
        Positions = new int[FS_PAIR];
        Orientations = new byte[FS_PAIR];
        Source.extractLocationData(Chromosomes, Positions, Orientations);

        Transcripts = new List[FS_PAIR];
        Transcripts[FS_UP] = Lists.newArrayList();
        Transcripts[FS_DOWN] = Lists.newArrayList();
        Source.extractTranscriptNames(Transcripts[FS_UP], Transcripts[FS_DOWN]);
    }

    public boolean isFusion() { return Source.VariantType.isFusion(); }
    public boolean isPointMutation() { return Source.VariantType.isPointMutation(); }

    public boolean singleGene() { return Source.GeneIds[FS_UP].equals(Source.GeneIds[FS_DOWN]); }

    public void setOrientation(final byte geneStrand)
    {
        if(geneStrand == POS_STRAND)
        {
            Orientations[FS_UP] = POS_ORIENT;
            Orientations[FS_DOWN] = NEG_ORIENT;
        }
        else
        {
            Orientations[FS_UP] = NEG_ORIENT;
            Orientations[FS_DOWN] = POS_ORIENT;
        }
    }

    public int[] getCodingBaseRange(int stream)
    {
        int codingBaseEnd = Source.CodingBasePositions[stream];
        int codingBaseLength = Source.CodingBases[stream].length();

        if((Orientations[stream] == POS_ORIENT) == (stream == FS_UP))
        {
            return new int[] { codingBaseEnd, codingBaseEnd + codingBaseLength - 1 };
        }
        else
        {
            return new int[] { codingBaseEnd - (codingBaseLength - 1), codingBaseEnd };
        }
    }

    public String getFullCodingBases(int streamPerspective, final int[] codingBaseRange)
    {
        if(!isFusion() || Orientations[FS_UP] != Orientations[FS_DOWN])
        {
            if(Orientations[FS_UP] == POS_ORIENT)
            {
                codingBaseRange[SE_START] = Source.CodingBasePositions[FS_UP];
                codingBaseRange[SE_END] = Source.CodingBasePositions[FS_DOWN];
                return Source.CodingBases[FS_UP] + Source.CodingBases[FS_DOWN];
            }
            else
            {
                codingBaseRange[SE_START] = Source.CodingBasePositions[FS_DOWN];
                codingBaseRange[SE_END] = Source.CodingBasePositions[FS_UP];
                return Source.CodingBases[FS_DOWN] + Source.CodingBases[FS_UP];
            }
        }

        int codingBaseEnd = Source.CodingBasePositions[streamPerspective];
        int codingBaseLength = Source.CodingBases[streamPerspective].length() - 1;

        if(streamPerspective == FS_UP)
        {
            if(Orientations[FS_UP] == POS_ORIENT)
            {
                codingBaseRange[SE_START] = codingBaseEnd;
                codingBaseRange[SE_END] = codingBaseEnd + codingBaseLength;
                return Source.CodingBases[FS_UP] + reverseStrandBases(Source.CodingBases[FS_DOWN]);
            }
            else
            {
                codingBaseRange[SE_START] = codingBaseEnd - codingBaseLength;
                codingBaseRange[SE_END] = codingBaseEnd;
                return reverseStrandBases(Source.CodingBases[FS_DOWN]) + Source.CodingBases[FS_UP];
            }
        }
        else
        {
            if(Orientations[FS_DOWN] == POS_ORIENT)
            {
                codingBaseRange[SE_START] = codingBaseEnd;
                codingBaseRange[SE_END] = codingBaseEnd + codingBaseLength;
                return Source.CodingBases[FS_DOWN] + reverseStrandBases(Source.CodingBases[FS_UP]);
            }
            else
            {
                codingBaseRange[SE_START] = codingBaseEnd - codingBaseLength;
                codingBaseRange[SE_END] = codingBaseEnd;
                return reverseStrandBases(Source.CodingBases[FS_UP]) + Source.CodingBases[FS_DOWN];
            }
        }
    }

}
