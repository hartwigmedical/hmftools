package com.hartwig.hmftools.sage.context;

import org.jetbrains.annotations.NotNull;

public class ReadContextDistance {

    private final String difference;
    private final int distance;

    public ReadContextDistance(final byte[] readBytes, final byte[] refBytes) {
        assert (readBytes.length == refBytes.length);

        int distance = 0;
        char[] diffChar = new char[readBytes.length];
        for (int i = 0; i < readBytes.length; i++) {
            byte refByte = refBytes[i];
            byte readByte = readBytes[i];

            if (refByte == readByte) {
                diffChar[i] = '.';
            } else {
                diffChar[i] = 'X';
                distance++;
            }

        }

        this.distance = distance;
        difference = new String(diffChar);
    }

    @NotNull
    public String difference() {
        return difference;
    }

    public int distance() {
        return distance;
    }
}
