package com.hartwig.hmftools.isofox.fusion;

import static com.hartwig.hmftools.common.bam.CigarUtils.cigarFromStr;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.N;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bam.ClippedSide;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.isofox.common.ReadRecord;

import htsjdk.samtools.Cigar;

public class SupplementaryJunctionData
{
    public final List<String> ReadIds;

    public int LocalJunctionPos;
    public byte LocalJunctionOrient;

    public String RemoteChromosome;
    public int RemoteJunctionPos;
    public byte RemoteJunctionOrient;

    public int MatchCount;

    public SupplementaryJunctionData(final String readId)
    {
        ReadIds = Lists.newArrayList(readId);

        LocalJunctionPos = 0;
        LocalJunctionOrient = 0;
        RemoteChromosome = "";
        RemoteJunctionPos = 0;
        RemoteJunctionOrient = 0;
        MatchCount = 0;
    }

    public boolean matches(final SupplementaryJunctionData other)
    {
        return LocalJunctionPos == other.LocalJunctionPos
                && RemoteChromosome.equals(other.RemoteChromosome)
                && RemoteJunctionPos == other.RemoteJunctionPos
                && LocalJunctionOrient == other.LocalJunctionOrient;
    }

    public static SupplementaryJunctionData fromReads(final ChimericReadGroup readGroup)
    {
        ReadRecord read = readGroup.reads().stream().filter(x -> x.hasSuppAlignment()).findFirst().orElse(null);
        if(read == null)
            return null;

        SupplementaryReadData suppData = SupplementaryReadData.extractAlignment(read.getSuppAlignment());

        if(suppData == null)
            return null;

        // find the junction from this read's SC and same for the supp mapping data
        ClippedSide scSide = ReadRecord.clippedSide(read);

        SupplementaryJunctionData suppJuncData = new SupplementaryJunctionData(read.Id);

        if(scSide.isLeft())
        {
            suppJuncData.LocalJunctionPos = read.getCoordsBoundary(SE_START);
            suppJuncData.LocalJunctionOrient = NEG_ORIENT;
        }
        else
        {
            suppJuncData.LocalJunctionPos = read.getCoordsBoundary(SE_END);
            suppJuncData.LocalJunctionOrient = POS_ORIENT;
        }

        suppJuncData.RemoteChromosome = suppData.Chromosome;
        Cigar remoteCigar = cigarFromStr(suppData.Cigar);
        scSide = ClippedSide.fromCigar(remoteCigar, true);

        if(scSide == null)
            return null;

        if(scSide.isLeft())
        {
            suppJuncData.RemoteJunctionPos = suppData.Position;
        }
        else
        {
            int skippedBases = remoteCigar.getCigarElements().stream()
                    .filter(x -> x.getOperator() == N || x.getOperator() == M || x.getOperator() == D)
                    .mapToInt(x -> x.getLength()).sum();

            suppJuncData.RemoteJunctionPos = suppData.Position + skippedBases - 1;
        }

        return suppJuncData;
    }

    public String toString()
    {
        return String.format("local(%d:%d) remote(%s:%d:%d) matched(%d)",
                LocalJunctionPos, LocalJunctionOrient, RemoteChromosome, RemoteJunctionPos, RemoteJunctionOrient, MatchCount);
    }

    public static void cacheSupplementaryJunctionCandidate(
            final ChimericReadGroup readGroup, final Map<Integer,List<SupplementaryJunctionData>> supplementaryJunctions)
    {
        SupplementaryJunctionData suppJuncData = fromReads(readGroup);

        if(suppJuncData == null)
            return;

        List<SupplementaryJunctionData> juncsByPos = supplementaryJunctions.get(suppJuncData.LocalJunctionPos);
        if(juncsByPos == null)
        {
            supplementaryJunctions.put(suppJuncData.LocalJunctionPos, Lists.newArrayList(suppJuncData));
            return;
        }

        SupplementaryJunctionData matchedData = juncsByPos.stream().filter(x -> x.matches(suppJuncData)).findFirst().orElse(null);

        if(matchedData != null)
        {
            ++matchedData.MatchCount;
        }
        else
        {
            juncsByPos.add(suppJuncData);
        }
    }

}
