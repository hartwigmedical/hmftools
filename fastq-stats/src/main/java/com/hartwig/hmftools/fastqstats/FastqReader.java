package com.hartwig.hmftools.fastqstats;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;

import org.jetbrains.annotations.NotNull;

class FastqReader {

    @NotNull
    private final BufferedInputStream reader;
    private final int size;

    FastqReader(@NotNull final InputStream inputStream) {
        this.size = 8192;
        this.reader = new BufferedInputStream(inputStream);
    }

    FastqReader(@NotNull final InputStream inputStream, final int size) {
        this.size = size;
        this.reader = new BufferedInputStream(inputStream, size);
    }

    @NotNull
    FastqData read() throws IOException {
        long yield = 0;
        long q30 = 0;
        int lineCount = 0;
        byte[] buf = new byte[size];
        int read;
        byte lastRead = 0;
        while ((read = reader.read(buf, 0, size)) != -1) {
            for (int i = 0; i < read; i++) {
                if (lastRead == '\r' && buf[i] == '\n') {
                    lastRead = buf[i];
                    continue;
                }
                if (buf[i] == '\r' || buf[i] == '\n') {
                    lastRead = buf[i];
                    lineCount++;
                    if (lineCount == 4) {
                        lineCount = 0;
                    }
                    continue;
                }
                if (lineCount == 3) {
                    yield++;
                    if (buf[i] >= 63) {
                        q30++;
                    }
                }
            }
        }
        return new FastqData(yield, q30);
    }

    void close() throws IOException {
        reader.close();
    }
}
