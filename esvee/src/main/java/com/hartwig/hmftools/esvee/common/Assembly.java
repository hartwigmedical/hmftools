package com.hartwig.hmftools.esvee.common;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

import com.hartwig.hmftools.esvee.html.DiagramSet;
import com.hartwig.hmftools.esvee.sequence.Sequence;
import com.hartwig.hmftools.esvee.sequence.SupportedAssembly;
import com.hartwig.hmftools.esvee.assembly.SequenceDecomposer;

import org.jetbrains.annotations.Nullable;

public abstract class Assembly implements Sequence
{
    public final String Name;

    public final String Assembly;
    public final byte[] AssemblyBases;

    protected final Map<Class<?>, List<Object>> mErrata = new HashMap<>();

    @Nullable
    private List<SequenceDecomposer.Node> mDecomposition = null;

    public Assembly(final String name, final String assembly)
    {
        Name = name;
        Assembly = assembly;
        AssemblyBases = assembly.getBytes();
    }

    @SuppressWarnings("unchecked")
    public <T> List<T> getAllErrata(final Class<T> type)
    {
        return Objects.requireNonNullElse((List<T>) mErrata.get(type), List.of());
    }

    public Map<Class<?>, List<Object>> getAllErrata()
    {
        return mErrata;
    }

    public void addErrata(final Object object)
    {
        mErrata.computeIfAbsent(object.getClass(), __ -> new ArrayList<>()).add(object);
    }

    public void addErrata(final List<Object> objects)
    {
        if(objects.isEmpty())
            return;

        mErrata.computeIfAbsent(objects.get(0).getClass(), __ -> new ArrayList<>()).addAll(objects);
    }

    public void addErrata(final Map<Class<?>, List<Object>> objects)
    {
        if(objects.isEmpty())
            return;

        objects.values().forEach(this::addErrata);
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

    public List<DiagramSet> getDiagrams()
    {
        return List.of();
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
