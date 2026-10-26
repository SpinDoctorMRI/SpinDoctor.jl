<!-- A card of the examples gallery: an image, the title and a description -->
<script setup lang="ts">
import { withBase } from 'vitepress'
export interface Props {
    href: string;
    src: string;
    caption: string;
    desc: string;
}

defineProps<Props>();
</script>

<template>
    <a class="card" :href="href">
        <img :src="withBase(src)" alt="">
        <div class="text">
            <h3>{{ caption }}</h3>
            <p>{{ desc }}</p>
        </div>
    </a>
</template>

<style scoped>
.card {
    display: flex;
    flex-direction: column;
    width: 250px;
    margin: 12px;
    border: 1px solid var(--vp-c-divider);
    border-radius: 12px;
    overflow: hidden;
    color: var(--vp-c-text-1);
    background-color: var(--vp-c-bg-soft);
    text-decoration: none;
    transition: border-color 0.2s ease, transform 0.2s ease;
}

.card:hover {
    border-color: var(--vp-c-brand-1);
    transform: translateY(-3px);
}

.card img {
    width: 100%;
    height: 180px;
    object-fit: contain;
    background-color: white;
}

.text {
    padding: 10px 14px 14px;
}

.text h3 {
    margin: 0 0 4px;
    padding: 0;
    border: 0;
    font-size: 1.05em;
}

.text p {
    margin: 0;
    font-size: 0.85em;
    line-height: 1.4;
    color: var(--vp-c-text-2);
}
</style>
