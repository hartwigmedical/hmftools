package com.hartwig.hmftools.esvee.common;

import java.util.List;

import com.hartwig.hmftools.esvee.old.Sequence;
import com.hartwig.hmftools.esvee.old.SupportedAssembly;
import com.hartwig.hmftools.esvee.old.SequenceDecomposer;

import org.jetbrains.annotations.Nullable;

public abstract class Assembly implements Sequence
{
    public final String Name;

    public final String Assembly;
    public final byte[] AssemblyBases;

    @Nullable
    private List<SequenceDecomposer.Node> mDecomposition = null;

    public Assembly(final String name, final String assembly)
    {
        Name = name;
        Assembly = assembly;
        AssemblyBases = assembly.getBytes();
    }

    @Override
    public String getName()
    {
        return Name;
    }

    @Override
    public String getBasesString()
    {
        return Assembly;
    }

    @Override
    public byte[] getBases()
    {
        return AssemblyBases;
    }

    @Override
    public abstract byte[] getBaseQuality();

    @Override
    public int getLength()
    {
        return Assembly.length();
    }

    @Override
    public List<SequenceDecomposer.Node> decompose()
    {
        if(mDecomposition == null)
        {
            if(this instanceof SupportedAssembly)
                mDecomposition = SequenceDecomposer.decompose((SupportedAssembly) this);
            else
                mDecomposition = SequenceDecomposer.decompose(getBases(), getBaseQuality());
        }
        return mDecomposition;
    }

    public void markDecompositionStale()
    {
        mDecomposition = null;
    }
}
