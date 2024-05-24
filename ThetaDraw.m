function ThetaDraw
    Source = load('genH2mat2_theta_first_test.dat','-ascii');
    theta_1 = Source(:,4);
    theta_2 = Source(:,5);
    tticks = Source(:,6);
    plot(tticks,theta_1,'b',tticks,theta_2,'r--');
    xlabel('t');
    ylabel('Theta_1,Theta_2');
end