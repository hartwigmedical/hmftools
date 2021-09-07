package com.hartwig.hmftools.linx.chaining;

import java.util.Collection;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.linx.types.SvVarData;

public class SvChainConnections
{
    private final Map<Integer,ChainState> mSvConnectionsMap;

    public SvChainConnections()
    {
        mSvConnectionsMap = Maps.newHashMap();
    }

    public void clear() { mSvConnectionsMap.clear(); }

    public void add(final SvVarData var, final ChainState state) { mSvConnectionsMap.put(var.id(), state); }

    public void remove(final SvVarData var) { mSvConnectionsMap.remove(var.id()); }

    public ChainState get(final SvVarData var) { return mSvConnectionsMap.get(var.id()); }

    public Collection<ChainState> values() { return mSvConnectionsMap.values(); }

    public int size() { return mSvConnectionsMap.size(); }

}
