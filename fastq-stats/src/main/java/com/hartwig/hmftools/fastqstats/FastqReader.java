package com.hartwig.hmftools.fastqstats;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;

public class FastqReader {
    private final BufferedInputStream reader;
    private final int size;

    public FastqReader(InputStream in){
        size = 8192;
        reader = new BufferedInputStream(in);

    }

    public FastqReader(InputStream in, int size) {
        this.size = size;
        reader = new BufferedInputStream(in, size);
    }

    /**
     * Reads Yield and Q30 counts from this input stream and returns them in a FastqData object
     * @return FastqData object containing yield and q30 counts
     * @throws IOException
     */
    public FastqData read() throws IOException{
        long yield = 0;
        long q30 = 0;
        int lineCount = 0;
        byte[] buf = new byte[size];
        int read;
        byte lastRead = 0;
        while((read = reader.read(buf, 0, size)) != -1){
            for(int i = 0; i < read; i ++){
                if(lastRead == '\r' && buf[i] == '\n'){
                    lastRead = buf[i];
                    continue;
                }
                if(buf[i] == '\r' || buf[i] == '\n'){
                    lastRead = buf[i];
                    lineCount ++;
                    if(lineCount == 4){
                        lineCount = 0;
                    }
                    continue;
                }
                if(lineCount == 3) {
                    yield ++;
                    if(buf[i] >= 63) {
                        q30 ++;
                    }
                }
            }
        }
        return new FastqData(yield, q30);
    }

    public void close() throws IOException {
        reader.close();
    }
}
