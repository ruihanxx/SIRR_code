close all
colors = ['#A2A9DD','#C4A5DE','#F0988C',"#EEB66D","#D94F33","#834026"];
colors = ["#DA5419",'#EDB120','#0072BD'];
colors={[0.85,0.33,0.10],[0.93,0.69,0.13],[0,0.45,0.74]};
close all
for idx2 = 1:length(res_sizes)
    res_size = res_sizes(idx2);
    for idx1 = length(conds):-1:1
        condA = -conds(idx1);
        condA = 10^condA;
        for j = 1:3
            figure(j+3*(idx2-1));
            %set(gca,'position',[0.1 0.1 0.8 0.92]);
            %axis normal;
            method1_idx = method1_result{idx1,idx2};
            method2_idx = method2_result{idx1,idx2};
            method3_idx = method3_result{idx1,idx2};
            qr_idx = qr_result{idx1,idx2};
            xlim([0 1+(J1-1)*J^(t2-1)*K])
            color = cell2mat(colors(idx1));
            a = semilogy(1+(0:J1-1)*J^(t2-1)*K,method1_idx(j,:),'*--', 'LineWidth', 2,...
                'Color', color); hold on
            a.Color(4)=0.5;
            semilogy((1:100)*5,method3_idx(j,:),'-', 'LineWidth', 2,...
                'Color', color); hold on
            sample = method2_idx(j,end-20:end);
            yline(mean(sample),'-', 'LineWidth', 2,'Color', color,'Alpha',0.5); hold on
            yline(qr_idx(j),':', 'LineWidth', 3, 'Color', color)
            % semilogy(method3_idx(j,:),'--', 'LineWidth', 2,...
            %     'Color', colors(idx1)); hold on
            patch([0, 1+(J1-1)*J^(t2-1)*K,1+(J1-1)*J^(t2-1)*K, 0], [quantile(sample,0.1), quantile(sample,0.1),quantile(sample,0.9), quantile(sample,0.9)],color, 'FaceAlpha', 0.1, 'EdgeColor', 'none')
            %legend(legend_list)
            if j == 1
                ylabel('Forward Error','FontSize',14);
                legend({'','',legend_cond{1},'','','','',legend_cond{2},'','','','',legend_cond{3},'',''},'Location','best')
            end
            if j == 2
                ylabel('Residual Error','FontSize',14);
            end
            if j == 3
                ylabel('Backward Error','FontSize',14);
            end
            xlabel('Iteration i');

        end
    end
end
%saveas(gca,'final.pdf')