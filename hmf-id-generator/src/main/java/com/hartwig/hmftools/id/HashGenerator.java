package com.hartwig.hmftools.id;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.bouncycastle.jcajce.provider.digest.SHA3;
import org.bouncycastle.util.encoders.Hex;

public class HashGenerator
{
    private final String mPassword;

    private static final List<String> PREFIXES = Lists.newArrayList("WIDE", "CPCT", "DRUP");
    private static final List<String> LOCATIONS = Lists.newArrayList("01", "02");
    private static final List<String> SUFFIXES = Lists.newArrayList("T", "TI", "TII", "TIII", "TIV");

    public HashGenerator(final String password)
    {
        mPassword = password;
    }

    public String hash(final String plaintext)
    {
        SHA3.Digest256 sha3 = new SHA3.Digest256();
        byte[] bytes = (plaintext + mPassword).getBytes();
        sha3.update(bytes);
        return Hex.toHexString(sha3.digest());
    }

    public static Map<String,String> precomputeHashes()
    {
        Map<String,String> hashMap = Maps.newHashMap();

        return null;
    }

    /*
    private fun precomputeHashes(generator: IdGenerator): Map<String, String> {
    val result = mutableMapOf<String, String>()
    val prefixes = setOf("WIDE", "CPCT", "DRUP")
    val locations = setOf("01", "02")
    val suffixes: Set<String> = setOf("T", "TI", "TII", "TIII", "TIV")

    for (prefix in prefixes) {
        logger.debug("  using prefix $prefix")
        for (location in locations) {
            logger.debug("   using location $location")
            for (suffix in suffixes) {
                logger.debug("     using suffix $suffix")
                for (i in 1..500000) {
                    val sample = prefix + location + i.toString().padStart(6, '0') + suffix
                    val newHash = generator.hash(sample)
                    result[newHash] = sample
                }
            }

        }
    }

    return result
}
     */

}
