package com.hartwig.hmftools.id;

import static java.lang.String.format;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.bouncycastle.jcajce.provider.digest.SHA3;
import org.bouncycastle.util.encoders.Hex;

public class HashGenerator
{
    private final String mPassword;
    private final int mMaxSampleCount;

    private static final List<String> PREFIXES = Lists.newArrayList("WIDE", "CPCT", "DRUP");
    private static final List<String> LOCATIONS = Lists.newArrayList("01", "02");
    private static final List<String> SUFFIXES = Lists.newArrayList("T", "TI", "TII", "TIII", "TIV");

    public HashGenerator(final String password, int maxSampleCount)
    {
        mPassword = password;
        mMaxSampleCount = maxSampleCount;
    }

    public Map<String,String> precomputeHashes()
    {
        // creates a pre-computed map of hash to expected original sample IDs
        // seems flawed since does not cover many of the new sample ID formats - does this matter?
        Map<String,String> hashMap = Maps.newHashMap();

        for(String prefix : PREFIXES)
        {
            for(String location : LOCATIONS)
            {
                for(String suffix : SUFFIXES)
                {
                    for(int i = 0; i < mMaxSampleCount; ++i)
                    {
                        String sample = format("%s%s%06d%s", prefix, location, i, suffix);
                        String newHash = hash(sample);
                        hashMap.put(newHash, sample);
                    }
                }
            }
        }

        return null;
    }

    public String hash(final String plaintext)
    {
        SHA3.Digest256 sha3 = new SHA3.Digest256();
        byte[] bytes = (plaintext + mPassword).getBytes();
        sha3.update(bytes);
        return Hex.toHexString(sha3.digest());
    }
}
