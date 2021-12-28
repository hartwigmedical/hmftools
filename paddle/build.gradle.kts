
plugins {
    id("com.hartwig.kotlin-conventions")
}

description = "HMF Tools - PADDLE"

dependencies {
    implementation(project(":hmf-common"))
    implementation(project(":patient-db"))
    implementation(libs.bcprov.jdk15on)
    testImplementation(libs.junit)
}
