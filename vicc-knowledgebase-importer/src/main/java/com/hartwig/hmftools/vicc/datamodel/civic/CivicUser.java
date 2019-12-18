package com.hartwig.hmftools.vicc.datamodel.civic;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicUser {

    @NotNull
    public abstract String username();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String displayName();

    @NotNull
    public abstract String role();

    @NotNull
    public abstract CivicOrganization organization();

    @Nullable
    public abstract String affiliation();

    @NotNull
    public abstract String featuredExpert();

    @Nullable
    public abstract String areaOfExpertise();

    @Nullable
    public abstract String bio();

    @Nullable
    public abstract String url();

    @NotNull
    public abstract String createdAt();

    @Nullable
    public abstract String lastSeenAt();

    @NotNull
    public abstract CivicAvatars avatars();

    @NotNull
    public abstract String avatarUrl();

    @Nullable
    public abstract String twitterHandle();

    @Nullable
    public abstract String facebookProfile();

    @Nullable
    public abstract String linkedinProfile();

    @Nullable
    public abstract String orcid();

    @Nullable
    public abstract String signupComplete();

    @Nullable
    public abstract String acceptedLicense();

    @NotNull
    public abstract String id();
}
