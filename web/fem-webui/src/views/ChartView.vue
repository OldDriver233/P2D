<template>
    <v-chart class="chart" :option="option" />
  </template>
  
<script setup lang="ts">
import { use } from "echarts/core";
import { CanvasRenderer } from "echarts/renderers";
import { LineChart } from "echarts/charts";
import {
  GridComponent
} from "echarts/components";
import VChart, { THEME_KEY } from "vue-echarts";
import { ref, provide } from "vue";
import axios from "axios";
import type { Ref } from "vue";
import type { ECBasicOption } from "echarts/types/dist/shared";
  
use([
  CanvasRenderer,
  LineChart,
  GridComponent
]);
  
provide(THEME_KEY, "light");

const x_axis: Ref<number[]> = ref([])
const data_points: Ref<number[]> = ref([])

axios.get("http://127.0.0.1:3000/api/result")
  .then(function (resp) {
    data_points.value = resp.data
    for(var i = 0; i < data_points.value.length; i++) {
      x_axis.value.push(i)
    }
    console.log(x_axis.value)
  })
  
const option: Ref<ECBasicOption> = ref({
  xAxis: {
    data: x_axis
  },
  yAxis: {
    type: 'value'
  },
  series: [
    {
      data: data_points,
      type: 'line',
      smooth: true,
      showSymbol: false
    }
  ],
  itemStyle: {

  }
});
</script>
  
  <style scoped>
  .chart {
    height: 400px;
  }
  </style>