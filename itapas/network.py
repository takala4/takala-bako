import numpy as np
import pandas as pd
from collections import defaultdict


class Network:
    """GMNS Plus フォーマットのネットワークデータを読み込み・管理するクラス"""

    def __init__(self, node_csv, link_csv, demand_csv):
        self._load_nodes(node_csv)
        self._load_links(link_csv)
        self._load_demand(demand_csv)
        self._build_adjacency()

    def _load_nodes(self, node_csv):
        df = pd.read_csv(node_csv)
        self.node_ids = df['node_id'].values
        self.num_nodes = len(self.node_ids)

        # ノードIDから内部インデックス (0-based) へのマッピング
        self.node_id_to_idx = {nid: i for i, nid in enumerate(self.node_ids)}
        self.idx_to_node_id = {i: nid for i, nid in enumerate(self.node_ids)}

        # ゾーンセントロイド (zone_id == node_id のノード)
        self.zone_ids = df.dropna(subset=['zone_id'])
        self.zone_ids = self.zone_ids[self.zone_ids['zone_id'] > 0]['zone_id'].astype(int).values

    def _load_links(self, link_csv):
        df = pd.read_csv(link_csv)
        self.num_links = len(df)

        self.link_from = np.array([self.node_id_to_idx[n] for n in df['from_node_id']])
        self.link_to = np.array([self.node_id_to_idx[n] for n in df['to_node_id']])

        # 容量: capacity (vehicles/hour) は1レーンあたりの値
        lanes = df['lanes'].values.astype(float)
        self.link_capacity = df['capacity'].values.astype(float) * lanes

        # 自由走行時間
        if 'vdf_fftt' in df.columns:
            self.link_fftt = df['vdf_fftt'].values.astype(float)
        else:
            # vdf_length_mi / vdf_free_speed_mph * 60 (分単位)
            self.link_fftt = df['vdf_length_mi'].values / df['vdf_free_speed_mph'].values * 60.0

        # BPRパラメータ
        self.link_alpha = df['vdf_alpha'].values.astype(float)
        self.link_beta = df['vdf_beta'].values.astype(float)

    def _load_demand(self, demand_csv):
        df = pd.read_csv(demand_csv)
        self.demand = {}
        self.origins = set()

        for _, row in df.iterrows():
            o = self.node_id_to_idx[int(row['o_zone_id'])]
            d = self.node_id_to_idx[int(row['d_zone_id'])]
            vol = float(row['volume'])
            if vol > 0:
                self.demand[(o, d)] = vol
                self.origins.add(o)

        self.origins = sorted(self.origins)

    def _build_adjacency(self):
        """隣接リスト (出リンク・入リンク) の構築"""
        self.out_links = defaultdict(list)
        self.in_links = defaultdict(list)

        for e in range(self.num_links):
            self.out_links[self.link_from[e]].append(e)
            self.in_links[self.link_to[e]].append(e)
