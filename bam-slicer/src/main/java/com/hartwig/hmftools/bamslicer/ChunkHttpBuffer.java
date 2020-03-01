package com.hartwig.hmftools.bamslicer;

import java.io.IOException;
import java.net.URL;
import java.util.List;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ExecutionException;

import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.util.concurrent.Futures;
import com.google.common.util.concurrent.ListenableFuture;
import com.google.common.util.concurrent.SettableFuture;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.Chunk;
import htsjdk.samtools.util.BlockCompressedFilePointerUtil;
import okhttp3.Call;
import okhttp3.Callback;
import okhttp3.Headers;
import okhttp3.OkHttpClient;
import okhttp3.Request;
import okhttp3.Response;
import okhttp3.ResponseBody;

class ChunkHttpBuffer {

    private static final Logger LOGGER = LogManager.getLogger(ChunkHttpBuffer.class);

    @NotNull
    private final LoadingCache<Long, ListenableFuture<byte[]>> chunkBuffer;
    private final int maxSize;
    @NotNull
    private final URL url;
    @NotNull
    private final ConcurrentSkipListMap<Long, Chunk> chunksPerOffset = new ConcurrentSkipListMap<>();
    @NotNull
    private final OkHttpClient httpClient;

    ChunkHttpBuffer(@NotNull OkHttpClient httpClient, @NotNull URL url, int maxSize, @NotNull List<Chunk> chunks) {
        this.httpClient = httpClient;
        this.url = url;
        this.maxSize = maxSize;

        for (Chunk chunk : chunks) {
            long chunkStart = BlockCompressedFilePointerUtil.getBlockAddress(chunk.getChunkStart());
            chunksPerOffset.put(chunkStart, chunk);
        }

        chunkBuffer = CacheBuilder.newBuilder().maximumSize(maxSize).build(new CacheLoader<Long, ListenableFuture<byte[]>>() {
            @Override
            @NotNull
            public ListenableFuture<byte[]> load(@NotNull Long offset) {
                Chunk chunkAtOffset = chunksPerOffset.get(offset);
                return getBytesForChunk(chunkAtOffset);
            }

            @Override
            public ListenableFuture<ListenableFuture<byte[]>> reload(@NotNull Long offset, @NotNull ListenableFuture<byte[]> oldBytes) {
                return Futures.immediateFuture(oldBytes);
            }
        });
    }

    @NotNull
    Pair<Long, byte[]> getEntryAtPosition(long position) throws IOException {
        long chunkOffset = chunksPerOffset.floorEntry(position).getKey();
        byte[] bytesAtOffset;
        try {
            if (chunkBuffer.getIfPresent(chunkOffset) == null) {
                refillBuffer(chunkOffset);
            }
            bytesAtOffset = chunkBuffer.getUnchecked(chunkOffset).get();
            chunkBuffer.invalidate(chunkOffset);
            if (chunkBuffer.size() < maxSize / 5 + 1) {
                refillBuffer(chunkOffset);
            }
            return Pair.of(chunkOffset, bytesAtOffset);
        } catch (InterruptedException | ExecutionException e) {
            throw new IOException("Could not read entry at position " + position + ". Cause: " + e.getMessage());
        }
    }

    private void refillBuffer(long chunkOffset) {
        int fillSize = (int) (maxSize * .75);
        chunksPerOffset.tailMap(chunkOffset, false).keySet().stream().limit(fillSize).forEach(key -> {
            if (chunkBuffer.getIfPresent(key) == null) {
                chunkBuffer.refresh(key);
            }
        });
    }

    @NotNull
    private ListenableFuture<byte[]> getBytesForChunk(@NotNull Chunk chunk) {
        long start = BlockCompressedFilePointerUtil.getBlockAddress(chunk.getChunkStart());
        long end = BlockCompressedFilePointerUtil.getBlockAddress(chunk.getChunkEnd());
        if (start <= end) {
            return readUrlBytes(start, end - start);
        } else {
            return Futures.immediateFailedFuture(new IllegalArgumentException("start offset is greater than end"));
        }
    }

    @NotNull
    private ListenableFuture<byte[]> readUrlBytes(long offset, long count) {
        Headers httpHeaders = new Headers.Builder().add("Range", "bytes=" + offset + "-" + (offset + count - 1)).build();
        Request request = new Request.Builder().url(url).headers(httpHeaders).build();
        SettableFuture<byte[]> bytesFuture = SettableFuture.create();
        httpClient.newCall(request).enqueue(retryingCallback(10, bytesFuture));
        return bytesFuture;
    }

    @NotNull
    private Callback retryingCallback(int retryCount, @NotNull SettableFuture<byte[]> resultFuture) {
        return new Callback() {
            @Override
            public void onFailure(@NotNull Call call, @NotNull IOException e) {
                retryCall(call, e, retryCount, resultFuture);
            }

            @Override
            public void onResponse(@NotNull Call call, @NotNull Response response) {
                ResponseBody body = response.body();
                try {
                    if (response.isSuccessful() && body != null) {
                        resultFuture.set(body.bytes());
                    } else {
                        String nullBody = body == null ? "body = null" : "";
                        Exception e = new IOException("Response " + response.code() + ": " + response.message() + "; " + nullBody);
                        retryCall(call, e, retryCount, resultFuture);
                    }
                } catch (final Exception e) {
                    retryCall(call, e, retryCount, resultFuture);
                } finally {
                    if (body != null) {
                        body.close();
                    }
                }
            }
        };
    }

    private void retryCall(@NotNull  Call call, @NotNull  Exception exception,  int remainingRetries,
            @NotNull  SettableFuture<byte[]> resultFuture) {
        if (remainingRetries <= 0) {
            LOGGER.error("Call {} [{}] failed with exception: {}",
                    call.request().method(),
                    call.request().headers(),
                    exception.getMessage());
            resultFuture.setException(exception);
        } else {
            call.clone().enqueue(retryingCallback(remainingRetries - 1, resultFuture));
        }
    }

    @NotNull
    URL url() {
        return url;
    }

    void closeHttpClient() {
        httpClient.dispatcher().executorService().shutdown();
    }
}
