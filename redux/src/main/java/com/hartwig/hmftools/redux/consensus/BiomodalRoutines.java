package com.hartwig.hmftools.redux.consensus;

import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;

public class BiomodalRoutines
{
    public static BaseQualPair determineBaseAndQual(
            final byte[] locationBases, final byte[] locationQuals, final String chromosome, int position, final RefGenome refGenome)
    {
        // TODO: use logic below to implement per base
        return null;
    }

    /*
    @Override
    public void buildConsensusRead(final List<SAMRecord> reads, final ConsensusState consensusState, final boolean hasIndels)
    {
        String chromosome = reads.get(0).getReferenceName();
        mChromosomeLength = mRefGenome.getChromosomeLength(chromosome);

        List<List<AnnotatedBase>> annotatedReads = annotateReads(reads);
        for(int i = 0; i < reads.size(); i++)
        {
            addModCAnnotation(reads.get(i), annotatedReads.get(i));
        }

        Collection<List<AnnotatedBase>> alignment = alignAnnotatedReads(annotatedReads);

        List<AnnotatedBase> consensusBases = getConsensusBases(
                alignment, records -> determineConsensus(chromosome, consensusState.IsForward, records));

        finalizeConsensusState(mRefGenome, reads, consensusState, hasIndels, consensusBases);

        SortedSet<Integer> modCReadIndices = Sets.newTreeSet();
        for(int i = 0; i < consensusBases.size(); i++)
        {
            if(consensusBases.get(i).Annotations.contains(MODC_ANNOTATION))
            {
                modCReadIndices.add(i);
            }
        }

        String mmValue = getMMValueFromModCReadIndices(consensusState.Bases, modCReadIndices, consensusState.IsForward);
        consensusState.Attributes.put(BASE_MODIFICATIONS_ATTRIBUTE, mmValue);
    }

    private AnnotatedBase determineConsensus(final String chromosome, boolean isForward, final List<AnnotatedBase> bases)
    {
        CigarOperator consensusOp = getConsensusCigarOp(bases);
        ExtendedRefPos pos = bases.get(0).Pos;
        byte[] locationBases = new byte[bases.size()];
        byte[] locationQuals = new byte[bases.size()];
        for(int i = 0; i < bases.size(); i++)
        {
            locationBases[i] = bases.get(i).Base;
            locationQuals[i] = bases.get(i).Qual;
        }

        int basePosition = pos.RefPos;
        if(consensusOp != M || basePosition < 1 || basePosition > mChromosomeLength)
        {
            basePosition = INVALID_POSITION;
        }

        BaseQualPair consensusBaseAndQual = mBaseBuilder.determineBaseAndQual(
                locationBases, locationQuals, chromosome, basePosition, false, null);

        byte consensusBase = !consensusBaseAndQual.isValid() ? ANY_BASE : consensusBaseAndQual.Base;
        byte consensusQual = BaseQualAdjustment.adjustBaseQual(consensusBaseAndQual.Qual);

        if(consensusQual == (byte) 0 && consensusOp == I)
        {
            return null;
        }

        AnnotatedBase consensus = new AnnotatedBase(pos, consensusBase, consensusQual, consensusOp);

        byte targetBase = isForward ? (byte) 'C' : (byte) swapDnaBase('C');
        if(consensusBase != targetBase)
        {
            return consensus;
        }

        int nonModCQualTotal = 0;
        int modCQualTotal = 0;
        for(AnnotatedBase annotatedBase : bases)
        {
            if(annotatedBase.Base != targetBase)
            {
                continue;
            }

            if(annotatedBase.Annotations.contains(MODC_ANNOTATION))
            {
                modCQualTotal += annotatedBase.Qual;
                continue;
            }

            nonModCQualTotal += annotatedBase.Qual;
        }

        if(modCQualTotal <= nonModCQualTotal)
        {
            return consensus;
        }

        consensus.Annotations.add(MODC_ANNOTATION);
        return consensus;
    }

    private static void addModCAnnotation(final SAMRecord read, List<AnnotatedBase> annotatedRead)
    {
        SortedSet<Integer> modCReadIndices = getModCReadIndices(read);
        for(int readIndex : modCReadIndices)
        {
            annotatedRead.get(readIndex).Annotations.add(MODC_ANNOTATION);
        }
    }
    */
}
